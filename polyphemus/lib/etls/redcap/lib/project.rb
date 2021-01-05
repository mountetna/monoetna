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
      models.each do |model_name, model_config|
        next unless build_model?(model_name.to_s)
        results[ model_name ] = {}

        results[ model_name ].update(
          Redcap::Model.create(
            model_name,
            model_config[:scripts],
            magma_models.model(model_name.to_s).template,
            template,
            dateshift_salt
          ).load(self)
        )
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

    def records(form, events)
      client.records(form, events)
    end

    private

    def build_model?(model_name)
      'all' == models_to_build ? true : models_to_build.include?(model_name)
    end

    def models_to_build
      config[:models_to_build]
    end

    def models
      config[:models]
    end

    def dateshift_salt
      config[:dateshift_salt]
    end

  end

end