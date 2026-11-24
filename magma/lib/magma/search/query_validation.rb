require_relative '../audit/project_to_template_mapping'

class Magma
  module Search
    # Rules Magma does not know: the saved-field shape, and a modality's template model.
    module QueryValidation
      COUNT_VERB = '::count'.freeze
      VERBS = { 'count' => [COUNT_VERB], 'filter' => ['::distinct', '::first'] }.freeze

      private

      def check_query!(query)
        return if query_array?(query) && query.first.is_a?(String)

        raise Etna::BadRequest, 'Query must be a non-empty JSON array whose first item is a model name, e.g. ["model_name", "::count"]'
      end

      def check_template_model!(name)
        return if template_model_names.include?(name)

        raise Etna::BadRequest, "Unknown model '#{name}'. Use a model from #{template_project}: #{template_model_names.sort.join(', ')}"
      end

      def result_verb(query)
        query.reverse.find { |part| verb?(part) }
      end

      # A json column hands back a JSONArray, which is not an ::Array.
      def query_array?(query)
        query.is_a?(Array) || query.is_a?(Sequel::Postgres::JSONArray)
      end

      def verb?(part)
        part.is_a?(String) && part.start_with?('::')
      end

      def template_model_names
        @template_model_names ||= Magma::ProjectToTemplateMapping.new(template_project).template_model_names
      end

      def template_project
        @template_project ||= Magma.instance.config(:template_project).to_s
      end
    end
  end
end
