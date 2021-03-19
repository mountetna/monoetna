require 'ostruct'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Magma
      class MaterializeDataWorkflow < Struct.new(
          :metis_client, :magma_client, :project_name,
          :model_name, :model_filters, :model_attributes_mask,
          :filesystem, :logger, :stub_files, :concurrency,
          :record_names, keyword_init: true)

        def initialize(**kwds)
          super(**({filesystem: Etna::Filesystem.new, concurrency: 10, record_names: "all"}.update(kwds)))
        end

        def magma_crud
          @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
        end

        def model_walker
          @model_walker ||= WalkModelTreeWorkflow.new(magma_crud: magma_crud, logger: logger)
        end

        def materialize_all(dest)
          templates = {}

          semaphore = Concurrent::Semaphore.new(concurrency)
          errors = Queue.new

          model_walker.walk_from(
              model_name,
              record_names,
              model_attributes_mask: model_attributes_mask,
              model_filters: model_filters,
              page_size: 20,
          ) do |template, document|
            logger&.info("Materializing #{template.name}##{document[template.identifier]}")
            templates[template.name] = template

            begin
              if (error = errors.pop(true))
                raise error
              end
            rescue ThreadError
            end

            semaphore.acquire
            Thread.new do
              begin
                materialize_record(dest, template, document)
              rescue => e
                errors << e
              ensure
                semaphore.release
              end
            end
          end

          semaphore.acquire(concurrency)

          begin
            if (error = errors.pop(true))
              raise error
            end
          rescue ThreadError
          end
        end

        def each_file(template, record, &block)
          results = []

          template.attributes.all.each do |attribute|
            if attribute.attribute_type == AttributeType::FILE_COLLECTION
              record[attribute.name]&.each_with_index do |file, i|
                results << [attribute, file, i]
              end
            elsif attribute.attribute_type == AttributeType::FILE
              results << [attribute, record[attribute.name], 0]
            end
          end

          results.each do |attr, file, idx|
            next if file.nil?
            next unless file.is_a?(Hash)
            next unless file['url']
            yield attr.name, file['url'], (file['original_filename'] || File.basename(file['path'])), idx
          end
        end

        def sync_metis_data_workflow
          @sync_metis_data_workflow ||= Etna::Clients::Metis::SyncMetisDataWorkflow.new(
              metis_client: metis_client,
              logger: logger,
              filesystem: filesystem)
        end

        def materialize_record(dest_dir, template, record)
          record_to_serialize = record.dup

          each_file(template, record) do |attr_name, url, filename, idx|
            if idx == 0
              record_to_serialize[attr_name] = []
            end

            dest_file = File.join(dest_dir, metadata_file_name(record_name: record[template.identifier], record_model_name: template.name, ext: "_#{attr_name}_#{idx}#{File.extname(filename)}"))
            sync_metis_data_workflow.copy_file(dest: dest_file, url: url, stub: stub_files)
            record_to_serialize[attr_name] << {file: dest_file, original_filename: filename}
          end

          dest_file = File.join(dest_dir, metadata_file_name(record_name: record[template.identifier], record_model_name: template.name, ext: '.json'))
          filesystem.mkdir_p(File.dirname(dest_file))
          json = record_to_serialize.to_json

          filesystem.with_writeable(dest_file, "w", size_hint: json.bytes.length) do |io|
            io.write(json)
          end
        end

        def metadata_file_name(record_name:, record_model_name:, ext:)
          "#{record_model_name}/#{record_name.gsub(/\s/, '_')}#{ext}"
        end
      end
    end
  end
end
