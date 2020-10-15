require 'ostruct'

module Etna
  module Clients
    class Magma
      class RecordSynchronizationWorkflow < Struct.new(:target_client, :source_client, :project_name, :ignore_update_errors, keyword_init: true)
        def target_models
          @target_models ||= begin
            target_client.retrieve(RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
          end
        end

        def crud
          @crud ||= MagmaCrudWorkflow.new(project_name: project_name, magma_client: target_client)
        end

        def source_models
          @source_models ||= source_client.retrieve(
              RetrievalRequest.new(
                  project_name: project_name, model_name: 'all')).models
        end

        # TODO: Add paging here to support large payloads
        def copy_model(model_name, copied_models = Set.new)
          if copied_models.include?(model_name)
            return
          end

          copied_models.add(model_name)

          descendants = source_models.to_directed_graph.descendants("project")
          (descendants[model_name] || []).each do |required_model|
            copy_model(required_model, copied_models) do |update|
              yield update if block_given?
            end
          end

          crud.page_records(model_name) do |documents|
            yield [model_name, documents] if block_given?
            begin
              crud.update_records do |update_request|
                documents.document_keys.each do |identifier|
                  record = documents.document(identifier)
                  record = record.each.select { |k, v| v != nil && !v.is_a?(Array) }.to_h
                  update_request.update_revision(model_name, identifier, record)
                end
              end
            rescue => e
              raise unless ignore_update_errors
            end
          end
        end

        def copy_all_models
          copied_models = Set.new
          target_models.model_keys.each do |model_name|
            copy_model(model_name, copied_models)
          end
        end
      end
    end
  end
end
