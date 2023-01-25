class Magma
  module WithFiltersModule
    def apply_filters(query)
      query = apply_filter_constraints(query)
      apply_filter_joins(query)
    end

    private

    def apply_filter_constraints(query)
      return query if @filters.empty?

      @filters.map do |filter|
        filter.flatten.map(&:constraint).concat(
          filter.flatten.map(&:subquery_constraints)).flatten
      end.inject(:+).each do |constraint|
        query = constraint.apply(query)
      end

      query
    end

    def apply_filter_joins(query)
      return query if @filters.empty?

      @filters.map do |filter|
        filter.flatten.map(&:join).flatten
      end.inject(:+).each do |join|
        query = join.apply(query)
      end

      query
    end
  end
end