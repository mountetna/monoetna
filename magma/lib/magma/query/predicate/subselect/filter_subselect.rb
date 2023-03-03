require_relative '../filter'

class Magma
  class FilterSubselectPredicate < Magma::FilterPredicate
    def self.verbs
      Magma::FilterPredicate.verbs.merge(@verbs)
    end

    verb do
      child do
        if Magma::SubqueryUtils.is_subquery?(self, @query_args)
          params = {
            predicate: self,
            question: @question,
            model_alias_name: @alias_name,
            query_args: @query_args
          }

          params[:subquery_class] = @parent_filter.subquery_config.magma_class if @parent_filter

          subquery = SubqueryFilter.new(**params)

          @subqueries << subquery

          subquery
        else
          RecordSubselectPredicate.new(@question, @model, @alias_name, *@query_args)
        end
      end
    end

    def create_filters
      invalid_argument!(@query_args.join(', ')) unless @query_args.all?{|q| q.is_a?(Array)}

      @filters = @query_args.map do |args|
        FilterSubselectPredicate.new(
          question: @question,
          model: @model,
          alias_name: @alias_name,
          parent_filter: self,
          query_args: args)
      end

      @query_args = []

      terminal TrueClass
    end
  end
end