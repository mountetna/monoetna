require_relative '../etls/redcap/redcap_etl_script_runner'

class Polyphemus
  class RedcapJob < Polyphemus::Job
    def self.as_json
      {
        name: "redcap",
        schema: Redcap::Loader.to_schema,
        params: {
          mode: [ 'default', 'strict', 'existing' ],
          model_names: 'string',
          commit: 'boolean'
        },
        secrets: [ :redcap_tokens ]
      }
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

    def run
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: request_params[:project_name],
        model_names: request_params[:model_names],
        redcap_tokens: request_params[:redcap_tokens],
        dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
        redcap_host: Polyphemus.instance.config(:redcap)[:host],
        magma_host: Polyphemus.instance.config(:magma)[:host],
        mode: request_params[:mode],
        config: request_params[:config]
      )

      magma_client = Etna::Clients::Magma.new(
        token: token,
        host: Polyphemus.instance.config(:magma)[:host])

      redcap_etl.run(magma_client: magma_client, commit: commit?, logger: $stdout)
    end

    private

    def array_or_string_param(param, allowed_values=["all"])
      request_params[param].is_a?(Array) || allowed_values.include?(request_params[param])
    end

    def commit?
      !!request_params[:commit]
    end
  end
end
