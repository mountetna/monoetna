module Redcap
  class Project
    attr_reader :config, :client, :magma_models, :magma_client, :logger, :project_name
    def initialize(token, project_name, config, magma_client, magma_models, logger)
      @client = Redcap::Client.new(config[:redcap_host], token)

      @config = config
      @project_name = project_name
      @magma_client = magma_client
      @magma_models = magma_models
      @logger = logger
    end

    def template
      @template ||= Redcap::Template.new(client.get_data(content: 'metadata'))
    end

    def fetch_records
      results = {}

      models.each do |model_name, model|
        results[ model_name ] = {}

        results[ model_name ].update(model.load)
      end

      results
    end

    def strict
      "strict" == config[:mode]
    end

    def eav_records
      @eav_records ||= client.get_record_eavs(field_opts(valid_fields)).reject do |record|
        record[:record] == 'test'
      end
    end

    def flat_records
      @flat_records ||= client.get_record_flat(
        field_opts([ 'record_id' ] + valid_fields)
      ).map do |record|
        record.merge(record: record[:record_id])
      end
    end

    private

    def valid_fields
      require 'pry'
      binding.pry
      return @valid_fields if @valid_fields

      all_fields = models.map do |model_name, model|
        model.scripts.map(&:fields).concat(model.identifier_fields)
      end.flatten.compact.uniq

      @valid_fields = all_fields & template.field_names
      invalid_fields = all_fields - template.field_names

      logger.write("The following fields are invalid: #{invalid_fields.join(", ")}") unless invalid_fields.empty?

      return @valid_fields
    end

    def field_opts(fields)
      fields.map.with_index do |f,i|
        ["fields[#{i}]", f ]
      end.to_h
    end

    def build_model?(model_name)
      'all' == models_to_build ? true : models_to_build.include?(model_name)
    end

    def models_to_build
      config[:models_to_build]
    end

    def etl_config
      config[:config].deep_symbolize_keys
    end

    def models
      @models ||= etl_config.map do |model_name, model_config|
        next unless build_model?(model_name.to_s)
        [
          model_name,
          Redcap::Model.create(
            self,
            model_name,
            model_config,
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
