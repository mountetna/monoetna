require_relative 'query_translator'
require_relative 'query_validation'
require_relative '../query/question'

class Magma
  class CrossProjectSearch
    include Magma::Search::QueryValidation

    MAX_QUERIES = 50

    def initialize(user)
      @user = user
    end

    def query(queries)
      check_queries!(queries)

      results = {}
      find_projects.each do |project_name|
        translator = Magma::Search::QueryTranslator.new(project_name)

        # One project's failure must not fail the batch, so errors become nil here.
        results[project_name] = queries.map do |query|
          translated = translator.translate(query)
          answer = launch_query(project_name, translated)
          clean_answer(query, answer)
        rescue Magma::Search::QueryTranslator::NotApplicable, Magma::QuestionError, ArgumentError, NameError
          nil
        rescue Sequel::DatabaseError => e
          Magma.instance.logger.log_error(e)
          nil
        end
      end

      { projects: results }
    end

    private

    def check_queries!(queries)
      raise Etna::BadRequest, 'queries must be an array' unless queries.is_a?(Array)
      raise Etna::BadRequest, 'queries must not be empty' if queries.empty?
      raise Etna::BadRequest, "at most #{MAX_QUERIES} queries per request" if queries.length > MAX_QUERIES

      queries.each_with_index do |query, index|
        check_query!(query)
      rescue Etna::BadRequest => e
        raise Etna::BadRequest, "Query #{index + 1}: #{e.message}"
      end
    end

    def find_projects
      Magma.instance.db[:models].distinct.select_map(:project_name).
        reject { |name| name.to_s.empty? || name == template_project }.
        sort
    end

    def launch_query(project_name, query)
      answer = Magma::Question.new(
        project_name,
        query,
        restrict: !@user.can_see_restricted?(project_name),
        user: @user,
        timeout: Magma.instance.config(:query_timeout)
      ).answer
      answer.respond_to?(:data) ? answer.data : answer
    end

    # Magma already answers in plain JSON (number, string, array); this only
    # drops blank values, trims strings and, for lists, de-duplicates.
    def clean_answer(query, value)
      verb = result_verb(query)
      return (value if value.is_a?(Numeric)) if verb == COUNT_VERB

      items = Array(value).map { |item| item.is_a?(String) ? item.strip : item }.reject { |item| item.nil? || item == '' }
      verb == '::first' ? items.first : items.uniq
    end
  end
end
