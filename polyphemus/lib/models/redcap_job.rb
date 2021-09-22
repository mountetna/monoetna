require "action_controller"

class Polyphemus
  class RedcapJob < Polyphemus::Job
    def self.as_json
      {
        name: "redcap",
        schema: Redcap::Loader.to_schema,
        secrets: [ :redcap_tokens ]
      }
    end

    def validate
      require_params(:redcap_tokens, :model_names)
      raise JobError, "redcap_tokens must be an array of tokens." unless request_params[:redcap_tokens].is_a?(Array)
      raise JobError, "model_names must be \"all\" or an array of model names." unless array_or_string_param(:model_names)

      if request_params[:mode]
        raise JobError, "mode must be nil, \"existing\", or \"strict\"." unless [nil, "existing", "strict"].include?(request_params[:mode])
      end
    end

    def run
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: request_params[:project_name],
        model_names: request_params[:model_names],
        redcap_tokens: request_params[:redcap_tokens],
        dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
        redcap_host: Polyphemus.instance.config(:redcap)[:host],
        magma_host: Polyphemus.instance.config(:magma)[:host],
        mode: request_params[:mode]
      )

      magma_client = Etna::Clients::Magma.new(
        token: user.token,
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
