require_relative '../etls/metis/loader'

class Polyphemus
  class MetisLinkerManifest < Polyphemus::WorkflowManifest
    def self.as_json
      {
        name: "metis",
        schema: Metis::Loader.to_schema,
        runtime_params: {
          commit: 'boolean'
        }
      }
    end

    def validate
    end

    private

    def array_or_string_param(param, allowed_values=["all"])
      request_params[param].is_a?(Array) || allowed_values.include?(request_params[param])
    end

  end
end
