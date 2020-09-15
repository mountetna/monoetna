require 'ostruct'

module Etna
  module Clients
    class Magma
      class MagmaCrudWorkflow < Struct.new(:magma_client, :project_name, :preview_only, keyword_init: true)
        attr_reader :recorded_updates

        def lookup_record(model_name, record_id)
          magma_client.retrieve(RetrievalRequest.new(project_name: project_name, record_names: [record_id], model_name: model_name))\
            .models.model(model_name).documents.document(record_id)
        end

        def update_records
          @recorded_updates ||= UpdateRequest.new(project_name: project_name)

          request = UpdateRequest.new(project_name: project_name)
          yield request
          @recorded_updates.revisions.update(request.revisions)

          magma_client.update(request) unless preview_only
        end
      end
    end
  end
end

