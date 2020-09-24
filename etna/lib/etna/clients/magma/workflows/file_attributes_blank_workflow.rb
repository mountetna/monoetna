require 'ostruct'
require 'pry'
require_relative './crud_workflow'

module Etna
  module Clients
    class Magma
      class FileAttributesBlankWorkflow < Struct.new(:magma_crud, :model_names, :project_name, keyword_init: true)

        def initialize(opts)
          super(**{}.update(opts))
        end

        def magma_client
          magma_crud.magma_client
        end

        def find_file_and_image_attributes
          {}.tap do |all_matches|
            model_names.each do |model_name|
              attribute_names = models.model(model_name).template.attributes.all.select {|a|
                a.attribute_type == Etna::Clients::Magma::AttributeType::FILE ||
                a.attribute_type == Etna::Clients::Magma::AttributeType::IMAGE}.map { |a|
                a.attribute_name}
              all_matches[model_name] = attribute_names
            end
          end
        end

        def each_revision
          find_file_and_image_attributes.each do |model_name, attribute_names|
            # Query all records in this model
            request = Etna::Clients::Magma::RetrievalRequest.new(project_name: project_name)
            request.model_name = model_name
            request.attribute_names = attribute_names
            request.record_names = 'all'
            documents = magma_client.retrieve(request).models.model(model_name).documents

            documents.document_keys.each do |record_name|
              document = documents.document(record_name)
              attribute_names.each do |attribute_name|
                yield [model_name, record_name, attribute_name] if !document[attribute_name]
              end
            end
          end
        end

        def set_file_attrs_blank
          magma_crud.update_records do |update_request|
            each_revision do |model_name, record_name, attribute_name|
              update_request.update_revision(model_name, record_name, revision_for(attribute_name))
            end
          end
        end

        def revision_for(attribute_name)
          {attribute_name => {path: "::blank"}}
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

