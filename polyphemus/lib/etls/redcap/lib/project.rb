module Redcap
  class Project
    attr_reader :config, :client, :magma_models, :logger
    def initialize(token, config, magma_models, logger)
      @client = Redcap::Client.new(config[:redcap_host], token)

      @config = config
      @magma_models = magma_models
      @logger = logger
    end

    def template
      @template ||= client.get_data(content: 'metadata')
    end

    def fetch_records
      results = {}

      models.each do |model_name, model|
        results[ model_name ] = {}

        results[ model_name ].update(model.load(self))
      end

      results
    end

    def forms
      @forms ||= template.map do |field|
        field[:form_name]
      end.uniq.map(&:to_sym)
    end

    def has_form?(form)
      forms.include?(form.to_sym)
    end

    def strict
      "strict" == config[:mode]
    end

    def eav_records
      @eav_records ||= client.get_record_eavs(
        fields: all_fields
      )
    end

    def flat_records
      @flat_records ||= client.get_record_flat(
        fields: [ 'record_id' ] + all_fields
      ).map do |record|
        [ record[:record_id], record ]
      end.to_h
    end

    def all_fields
      models.map do |model_name, model|
        model.scripts.map do |script|
          script.forms.values.map do |form|
            form.fields
          end
        end
      end.flatten
    end

    def records(fields, events, form_name, labels)
      eavs = eav_records.select do | eav|
        fields.include?(eav[:field_name])
      end

      return nil unless eavs

      records = eavs.group_by do |eav|
        events ?  [ eav[:record], eav[:redcap_event_name] ] : eav[:record]
      end.map do |record_id, record_eavs|
        [
          record_id,
          Redcap::Client::Record.new(
            record_eavs,
            form_name,
            flat_records[record_id],
            labels
          ).record
        ]
      end.to_h.compact

      return records
    end

    private

    def build_model?(model_name)
      'all' == models_to_build ? true : models_to_build.include?(model_name)
    end

    def models_to_build
      config[:models_to_build]
    end

    def models
      @models ||= config[:models].map do |model_name, model_config|
        next unless build_model?(model_name.to_s)
        [
          model_name,
          Redcap::Model.create(
            model_name,
            model_config[:scripts],
            magma_models.model(model_name.to_s).template,
            template,
            dateshift_salt
          )
        ]
      end.compact.to_h
    end

    def dateshift_salt
      config[:dateshift_salt]
    end
  end
end
