require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Magma
      class MaterializeDataWorkflow < Struct.new(
          :metis_client, :magma_client, :project_name,
          :model_name, :model_filters, :model_attributes_mask,
          :filesystem, :logger, :stub_files,
          :skip_tmpdir, keyword_init: true)

        def initialize(**kwds)
          super(**({filesystem: Etna::Filesystem.new}.update(kwds)))
        end

        def magma_crud
          @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
        end

        def model_walker
          @model_walker ||= WalkModelTreeWorkflow.new(magma_crud: magma_crud, logger: logger)
        end

        def materialize_all(dest)
          tmpdir = skip_tmpdir ? nil : filesystem.tmpdir

          begin
            model_walker.walk_from(
                model_name,
                model_attributes_mask: model_attributes_mask,
                model_filters: model_filters,
            ) do |template, document|
              logger&.info("Materializing #{template.name}##{document[template.identifier]}")
              materialize_record(dest, tmpdir, template, document)
            end
          ensure
            filesystem.rm_rf(tmpdir) unless skip_tmpdir
          end
        end

        def each_root_record
          request = RetrievalRequest.new(project_name: project_name, model_name: model_name, record_names: "all",
              filter: filter, page_size: 100, page: 1)
          magma_crud.page_records(model_name, request) do |response|
            model = response.models.model(model_name)
            template = model.template
            model.documents.document_keys.each do |key|
              yield template, model.documents.document(key)
            end
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
              skip_tmpdir: skip_tmpdir,
              filesystem: filesystem)
        end

        def materialize_record(dest_dir, tmpdir, template, record)
          record_to_serialize = record.dup

          each_file(template, record) do |attr_name, url, filename, idx|
            if idx == 0
              record_to_serialize[attr_name] = []
            end

            dest_file = File.join(dest_dir, metadata_file_name(record_name: record[template.identifier], record_model_name: template.name, ext: "_#{attr_name}_#{idx}#{File.extname(filename)}"))
            sync_metis_data_workflow.copy_file(bin_root_dir: dest_dir, tmpdir: tmpdir, dest: dest_file, url: url, stub: stub_files)
            record_to_serialize[attr_name] << { file: dest_file, original_filename: filename }
          end

          dest_file = File.join(dest_dir, metadata_file_name(record_name: record[template.identifier], record_model_name: template.name, ext: '.json'))
          filesystem.mkdir_p(File.dirname(dest_file))
          filesystem.with_writeable(dest_file, "w") do |io|
            io.write(record_to_serialize.to_json)
          end
        end

        def metadata_file_name(record_name:, record_model_name:, ext:)
          "#{record_model_name}/#{record_name.gsub(/\s/, '_')}#{ext}"
        end
      end
    end
  end
end
