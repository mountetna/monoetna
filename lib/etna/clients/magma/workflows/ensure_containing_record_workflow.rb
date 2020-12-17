require 'ostruct'
require_relative './crud_workflow'

module Etna
  module Clients
    class Magma
      class EnsureContainingRecordWorkflow < Struct.new(:magma_crud, :models,  keyword_init: true)
        def self.from_remote_models(magma_crud:, **opts)
          opts[:models] = magma_crud.magma_client.retrieve(RetrievalRequest.new(project_name: magma_crud.project_name, model_name: 'all')).models
          opts[:magma_crud] = magma_crud
          self.class.new(**opts)
        end

        def magma_client
          magma_crud.magma_client
        end

        def ensure_record(model_name, record_identifiers)
          raise "Could not find containing model #{model_name} defined." unless (model = models.model(model_name))

          raise "Identifiers #{record_identifiers} do not contain #{model_name}" unless (id = record_identifiers[model_name])
          record = magma_crud.lookup_record(model_name, id)

          if record.nil?
            attrs = { model.template.identifier => id }

            parent_attr = model.template.attributes.all.select { |a| a.attribute_type == AttributeType::PARENT }.first
            unless parent_attr.nil?
              parent_attribute_name = parent_attr.attribute_name
              parent_identifier = ensure_record(model.template.parent, record_identifiers)
              attrs.update({parent_attribute_name => parent_identifier})
            end

            magma_crud.update_records do |update_request|
              update_request.update_revision(model_name, id, attrs)
            end
          end

          id
        end
      end
    end
  end
end
