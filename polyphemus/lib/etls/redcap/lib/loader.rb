module Redcap
  class Loader
    attr_reader :records, :magma_models, :config, :logger
    def initialize(config, magma_client, logger=STDOUT)
      @config = config
      @records = {}
      @logger = logger

      # Fetch existing records separately, since tables
      #   don't have identifiers and including
      #   "record_names": "all", "attribute_names": "identifiers"
      #   will not return their
      #   templates.
      @magma_models = magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(
          project_name: @config[:project_name]
        )
      ).models
    end

    def run
      update_records_from_project

      patch_tables

      records
    end

    private

    def update_records_from_project
      tokens.each do |token|
        project = Redcap::Project.new(token, config, magma_models, logger)

        project.fetch_records.each do |model_name, model_records|
          records[model_name] = {} unless records.keys.include?(model_name)
          records[model_name].update(
            filter_records(model_name, model_records)
          )
        end
      end
    end

    def patch_tables
      # find any table attributes
      magma_models.model_keys.each do |model_name|
        magma_models.model(model_name.to_s).template.attributes.tap do |atts|
          atts.attribute_keys.each do |att_name|
            if atts.attribute(att_name).attribute_type == 'table'
              # look for corresponding entries in the database
              link_model_name = atts.attribute(att_name).link_model_name.to_sym

              if records[link_model_name]
                records[link_model_name].to_a.group_by do |(record_name,record)|
                  record[ model_name.to_sym ]
                end.each do |model_record_name, revisions|
                  temp_id_names = revisions.map(&:first)
                  records[ model_name.to_sym ] ||= {}
                  records[ model_name.to_sym ][ model_record_name ] ||= {}
                  records[ model_name.to_sym ][ model_record_name ][ att_name ] = temp_id_names
                end
              end
            end
          end
        end
      end
    end

    def tokens
      config[:tokens]
    end

    def filter_records(model_name, model_records)
      return model_records
      # return model_records unless restrict_records_to_update?

      model_records = filter_records_by_allowed_list(model_name, model_records)

      model_records
    end

    def filter_records_by_allowed_list(model_name, model_records)
      model_records.slice(*allowed_record_names(model_name))
    end

    def restrict_records_to_update?
      !!config[:records_to_update]
    end

    def allowed_record_names(model_name)
      return config[:records_to_update] unless "existing" == config[:records_to_update]

      # For tables, we actually need to filter on their parent
      #   model record_names....so this is broken for now.
      magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(
          project_name: @config[:project_name],
          model_name: model_name.to_s,
          attribute_names: "identifier",
          record_names: "all"
        )
      ).models.model(model_name.to_s).documents.document_keys
    end
  end
end
