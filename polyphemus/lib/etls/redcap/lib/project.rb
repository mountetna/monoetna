module Redcap
  class Project
    attr_reader :config, :client, :magma_models
    def initialize(client, config, magma_models)
      @client = client
      @config = config
      @magma_models = magma_models
    end

    def template
      @template ||= client.get_data(content: 'metadata')
    end

    def fetch_records
      results = {}
      models.each do |model_name, model_config|
        next if model_config[:disabled]
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

    def models
      config[:models]
    end

    def dateshift_salt
      config[:dateshift_salt]
    end

  end

end