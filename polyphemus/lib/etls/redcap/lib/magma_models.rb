module Redcap
  class MagmaModelsWrapper
    attr_reader :project_name, :magma_client
    def initialize(magma_client, project_name)
      @project_name = project_name
      @magma_client = magma_client

      # Fetch existing records separately, since tables
      #   don't have identifiers and including
      #   "record_names": "all", "attribute_names": "identifiers"
      #   will not return their
      #   templates.
      @magma_models = magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(
          project_name: project_name
        )
      ).models
    end

    def models
      @magma_models
    end

    def record_ids_for_model(model_name, model_records)
      # Fetches the relevant record ids for the given model.
      # For a non-table model, just requests the records directly.
      # For a table model, requests the parent's record names and then
      #   returns all the temporary ids from model_records
      #   that link to the parent records.
      return model_record_names(model_name) unless is_table?(model_name)

      parent_record_names = model_record_names(parent_model_name(model_name))

      model_records.select do |temp_id, record|
        parent_record_names.include?(record[parent_model_name(model_name).to_sym])
      end.keys
    end

    def parent_model_name(model_name)
      models.model(model_name.to_s).template.parent
    end

    def is_table?(model_name)
      # From a set of magma models, can only tell if a model
      #    is a table from its parent's definition of the
      #    child attribute.
      models.model(
        parent_model_name(model_name)
      ).template.attributes.attribute(
        model_name.to_s
      ).attribute_type == Etna::Clients::Magma::AttributeType::TABLE
    end

    def model_record_names(model_name)
      @record_names ||= {}
      @record_names[model_name.to_s] ||= magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(
          project_name: project_name,
          model_name: model_name.to_s,
          attribute_names: "identifier",
          record_names: "all"
        )
      ).models.model(model_name.to_s).documents.document_keys
    end
  end
end
