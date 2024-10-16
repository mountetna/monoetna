class Magma
  class FilterPredicate < Magma::Predicate
    attr_reader :model, :arguments

    def initialize question:, model:, alias_name:, query_args:, parent_filter: nil
      super(question)

      @filters = []
      @model = model
      @alias_name = alias_name
      @parent_filter = parent_filter
      process_args(query_args)
    end

    def create_filters
      invalid_argument!(@query_args.join(', ')) unless @query_args.all?{|q| q.is_a?(Array)}

      @filters = @query_args.map do |args|
        FilterPredicate.new(
          question: @question,
          model: @model,
          alias_name: @alias_name,
          parent_filter: self,
          query_args: args)
      end

      @query_args = []

      terminal TrueClass
    end

    verb '::or' do
      child :create_filters

      join :join_filters

      subquery_config Magma::SubqueryConfig.new(magma_class: Magma::SubqueryOuter)

      constraint do
        or_constraint(all_filter_constraints)
      end
    end

    verb '::and' do
      child :create_filters

      join :join_filters

      subquery_config Magma::SubqueryConfig.new(magma_class: Magma::SubqueryInner)

      constraint do
        and_constraint(all_filter_constraints)
      end
    end

    verb '::any' do
      subquery :join_subqueries

      subquery_config Magma::SubqueryConfig.new(condition: "> 0")
    end

    verb '::every' do
      subquery :join_subqueries

      subquery_config Magma::SubqueryConfig.new(condition: "= count(*)")
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
          RecordPredicate.new(@question, @model, @alias_name, *@query_args)
        end
      end
    end

    def subquery
      join_subqueries.concat(join_filter_subqueries)
    end

    def subquery_class
      @verb.do(:subquery_class)
    end

    def all_filter_constraints
      @filters.map do |filter|
        filter.flatten.map(&:constraint).concat(
          filter.flatten.map(&:subquery_constraints)).flatten
      end.concat(subquery_constraints).flatten.compact
    end
  end
end
