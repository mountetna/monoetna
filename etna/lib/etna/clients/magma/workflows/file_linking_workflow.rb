require 'ostruct'
require_relative './crud_workflow'

module Etna
  module Clients
    class Magma
      class FileLinkingWorkflow < Struct.new(:magma_crud, :model_name, :metis_client, :project_name, :bucket_name, :matching_expressions, :attribute_options, keyword_init: true)
        PATIENT_TIMEPOINT_REGEX = /([^-]+-[^-]+)-(DN?[0-9]+).*/

        def initialize(opts)
          super(**{attribute_options: {}, matching_expressions: []}.update(opts))
        end

        def magma_client
          magma_crud.magma_client
        end

        def find_matches
          {}.tap do |all_matches|
            metis_client.folders(
                project_name: project_name,
                bucket_name: bucket_name
            ).each do |folder|
              metis_client.list_folder(
                  Etna::Clients::Metis::ListFolderRequest.new(
                      project_name: project_name,
                      bucket_name: bucket_name,
                      folder_path: folder.folder_path,
                  ),
              ).files.all.each do |file|
                matches = matching_expressions
                          .map { |regex, attribute_name| [regex.match(file.file_path), regex, attribute_name] }
                          .select { |match, regex, attribute_name| !match.nil? }

                if matches.length > 1
                  raise "File path #{file.file_path} matches multiple regex, #{matches.map(&:second)}.  Please modify the matching expressions to disambiguate"
                end

                if matches.length == 1
                  match, _, attribute_name = matches.first
                  match_map = match.names.zip(match.captures).to_h
                  key = [match_map, attribute_name]

                  if attribute_options.dig(attribute_name, :file_collection)
                    (all_matches[key] ||= []).push(file.file_path)
                  else
                    if all_matches.include?(key)
                      raise "Field #{attribute_name} for #{match_map} identified for two files, #{file.file_path} and #{all_matches[key]}.  Please modify the existing matching expressionts to disambiguate"
                    end

                    all_matches[key] = [file.file_path]
                  end
                end
              end
            end
          end
        end

        # Subclasses should override this to implement custom logic for how regex matches should match to linking.
        def matches_to_record_identifiers(match_data)
          {"project" => project_name}
        end

        def patient_timepoint_from(word)
          match = PATIENT_TIMEPOINT_REGEX.match(word)
          return {} unless match

          return {
              'patient' => match[1],
              'timepoint' => "#{match[1]}-#{match[2]}",
          }
        end

        def each_revision
          find_matches.each do |key, file_paths|
            match_map, attribute_name = key
            record_identifiers = matches_to_record_identifiers(match_map)
            id = containing_record_workflow.ensure_record(model_name, record_identifiers)
            file_paths.each do |file_path|
              yield [id, revision_for(id, attribute_name, file_path, match_map, record_identifiers)]
            end
          end
        end

        def link_files
          magma_crud.update_records do |update_request|
            each_revision do |id, revision|
              update_request.update_revision(model_name, id, revision)
            end
          end
        end

        def revision_for(id, attribute_name, file_path, match_map, record_identifiers)
          if attribute_options.dig(attribute_name, :file_collection)
            file_path = ::File.dirname(file_path)
            {attribute_name => "https://metis.ucsf.edu/#{project_name}/browse/#{bucket_name}/#{file_path}"}
          else
            {attribute_name => {
              path: "metis://#{project_name}/#{bucket_name}/#{file_path}",
              original_filename: File.basename(file_path)}}
          end
        end

        def containing_record_workflow
          @containing_record_workflow ||= EnsureContainingRecordWorkflow.new(magma_crud: magma_crud, models: models)
        end

        def models
          @models ||= begin
            magma_client.retrieve(RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
          end
        end
      end
    end
  end
end

