require_relative '../etls/metis/loader'

class Polyphemus
  class MetisJob < Polyphemus::Job
    def self.as_json
      {
        name: "metis",
        schema: Metis::Loader.to_schema,
        params: { }
      }
    end

    def validate
    end

    def run
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
