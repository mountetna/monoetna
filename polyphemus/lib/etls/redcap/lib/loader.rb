module Redcap
  class Loader
    attr_reader :records, :magma_models_wrapper, :config, :logger
    def initialize(config, magma_models_wrapper, logger=STDOUT)
      @config = config
      @records = {}
      @logger = logger

      @magma_models_wrapper = magma_models_wrapper
    end

    def run
      update_records_from_project

      patch_tables

      return records, @records_to_blank
    end

    private

    def update_records_from_project
      tokens.each do |token|
        project = Redcap::Project.new(token, config, magma_models_wrapper.models, logger)

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
                if strict_update
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
      return model_records unless restrict_records_to_update?

      model_records = filter_records_by_allowed_list(model_name, model_records)

      model_records
    end

    def filter_records_by_allowed_list(model_name, model_records)
      model_records.slice(*allowed_record_names(model_name))
    end

    def restrict_records_to_update?
      !!config[:records_to_update] && !strict_update
    end

    def strict_update
      "strict" == config[:records_to_update]
    end

    def existing_records_update
      "existing" == config[:records_to_update]
    end

    def allowed_record_names(model_name)
      return config[:records_to_update] unless existing_records_update

      magma_models_wrapper.record_ids_for_model(model_name)
    end
  end
end
