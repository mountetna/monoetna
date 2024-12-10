class Polyphemus
  class RedcapLoaderManifest < Polyphemus::WorkflowManifest

    class << self
      def as_json
        {
          name: 'redcap',
          schema: Redcap::Loader.to_schema,
          runtime_params: {
            mode: [ {
              value: 'default',
              default: true,
              description: 'update existing and append new records'
            }, {
              value: 'strict',
              description: 'discard all existing records and append new records'
            }, {
              value: 'existing',
              description: 'only update existing records'
            } ],
            model_names: {
              type: 'options',
              value: 'model_names'
            },
            commit: 'boolean'
          },
          secrets: [ :redcap_tokens ]
        }
      end
    end

    def validate
      require_params(:redcap_tokens, :model_names)
      raise JobError, "redcap_tokens must be a comma-separated list of tokens." unless comma_separated?(:redcap_tokens, '[0-9A-F]{32}')
      raise JobError, "model_names must be \"all\" or a comma-separated list of model names." unless comma_separated?(:model_names, '[a-z]*(_[a-z]+)*')

      raise JobError, "mode must be \"default\", \"existing\", or \"strict\"." unless ["default", "existing", "strict"].include?(request_params[:mode])
    end

    def comma_separated?(param_name, pattern)
      request_params[param_name] =~ /\A#{pattern}(,\s*#{pattern})*\z/
    end

    private

    def array_or_string_param(param, allowed_values=["all"])
      request_params[param].is_a?(Array) || allowed_values.include?(request_params[param])
    end

  end
end
