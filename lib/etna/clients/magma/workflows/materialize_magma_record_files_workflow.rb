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
          keyword_init: true)

        def initialize(**kwds)
          super(**({filesystem: Etna::Filesystem.new}.update(kwds)))
        end

        def magma_crud
          @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
        end

        def model_walker
          @model_walker ||= WalkModelTreeWorkflow.new(magma_crud: magma_crud, logger: logger)
        end

        def with_materialized_dir(
            dir = filesystem.tmpdir,
            remove_on_failure: true,
            remove_on_success: true,
            &block)
          begin
            model_walker.walk_from(
                model_name,
                model_attributes_mask: model_attributes_mask,
                model_filters: model_filters,
            ) do |template, document|
              logger&.info("Materializing #{template.name}##{document[template.identifier]}")
              materialize_record(dir, template, document)
            end

            yield dir
            filesystem.rm_rf(dir) if remove_on_success
          rescue => e
            filesystem.rm_rf(dir) if remove_on_failure
            raise e
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

        def materialize_record(dest_dir, template, record)
          record_to_serialize = record.dup
          metadata_path = metadata_file_name(record_name: record[template.identifier], record_model_name: template.name)

          each_file(template, record) do |attr_name, url, filename, idx|
            metadata = metis_client.file_metadata(url)
            etag = metadata[:etag]
            size = metadata[:size]

            if idx == 0
              record_to_serialize[attr_name] = []
            end

            dest_file = bin_file_name(etag: etag)
            record_to_serialize[attr_name] << { file: dest_file, original_filename: filename }

            # Already materialized, continue
            if filesystem.exist?(dest_file)
              next
            end

            logger&.info("materializing file #{filename} (#{size} bytes)")
            filesystem.mkdir_p(File.dirname(File.join(dest_dir, dest_file)))

            upload_timings = []
            upload_amount = 0
            last_rate = 0.00001

            filesystem.with_writeable(File.join(dest_dir, dest_file), "w", size_hint: size) do |io|
              if stub_files
                io.write("(stub) #{filename}: #{size} bytes")
              else
                metis_client.download_file(url) do |chunk|
                  io.write(chunk)

                  upload_timings << [chunk.length, Time.now.to_f]
                  upload_amount += chunk.length

                  if upload_timings.length > 50
                    s, _ = upload_timings.shift
                    upload_amount -= s
                  end

                  _, start_time = upload_timings.first
                  _, end_time = upload_timings.last

                  if start_time == end_time
                    next
                  end

                  rate = upload_amount / (end_time - start_time)

                  if rate / last_rate > 1.2 || rate / last_rate < 0.8
                    logger&.info("Uploading #{Etna::Formatting.as_size(rate)} per second")

                    if rate == 0
                      last_rate = 0.0001
                    else
                      last_rate = rate
                    end
                  end
                end
              end
            end
          end

          dest_file = File.join(dest_dir, metadata_path)
          filesystem.mkdir_p(File.dirname(dest_file))
          filesystem.with_writeable(dest_file, "w") do |io|
            io.write(record_to_serialize.to_json)
          end
        end

        def metadata_file_name(record_name:, record_model_name:)
          "#{record_model_name}/#{record_name.gsub(/\s/, '_')}.json"
        end

        def bin_file_name(etag:)
          "bin/#{etag}"
        end
      end
    end
  end
end
