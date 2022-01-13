module Redcap
  class Loader
    def self.to_schema
      {
        "$schema": "http://json-schema.org/draft-07/schema#",
	"$id": "https://example.com/product.schema.json",
	title: "Redcap Loader",
	description: "This loader takes data from a Redcap host and imports it into Magma according to a specified mapping template. It requires a Redcap API token.",
        definitions: [
          Redcap::Model,
          Redcap::Script,
          Redcap::Entity,
          Redcap::Value
        ].map(&:to_schema).reduce(&:merge),
	type: "object",
        additionalProperties: {
          "$ref": "#/definitions/redcap_model"
        }
      }
    end

    attr_reader :records, :magma_models_wrapper, :config, :logger
    def initialize(config, project_name, magma_client, logger=STDOUT)
      @config = config
      @records = {}
      @logger = logger

      @project_name = project_name
      @magma_client = magma_client
      @magma_models_wrapper = Redcap::MagmaModelsWrapper.new(magma_client, project_name)
    end

    def run
      update_records_from_project

      patch_tables

      return records, @records_to_blank
    end

    private

    def update_records_from_project
      tokens.each do |token|
        project = Redcap::Project.new(
          token, @project_name, config, @magma_client, magma_models_wrapper.models, logger
        )
        
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
      magma_models_wrapper.models.model_keys.each do |model_name|
        magma_models_wrapper.models.model(model_name.to_s).template.attributes.tap do |atts|
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

                # Blank out records that have no data if user specifies a
                #   "strict" update
                if strict_mode
                  found_record_names = records[ model_name.to_sym ].keys.map(&:to_s)
                  @records_to_blank ||= {}
                  @records_to_blank[ model_name.to_sym ] ||= []
                  @records_to_blank[ model_name.to_sym ] = magma_models_wrapper.model_record_names(model_name).reject do |r|
                    found_record_names.include?(r)
                  end
                  @records_to_blank[ model_name.to_sym ].each do |model_record_name|
                    records[ model_name.to_sym ] ||= {}
                    records[ model_name.to_sym ][ model_record_name ] ||= {}
                    records[ model_name.to_sym ][ model_record_name ][ att_name ] = [] # empty array in Magma removes existing table records for an attribute
                  end
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
      return model_records unless restrict_mode?

      filter_records_by_allowed_list(model_name, model_records)
    end

    def filter_records_by_allowed_list(model_name, model_records)
      model_records.slice(*allowed_record_names(model_name, model_records))
    end

    def restrict_mode?
      !default_mode && !strict_mode && existing_mode
    end

    def default_mode
      "default" == config[:mode]
    end

    def strict_mode
      "strict" == config[:mode]
    end

    def existing_mode
      "existing" == config[:mode]
    end

    def allowed_record_names(model_name, model_records)
      magma_models_wrapper.record_ids_for_model(model_name, model_records)
    end
  end
end
