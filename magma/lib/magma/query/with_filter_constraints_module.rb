class Magma
  module WithFilterConstraintsModule
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
  end
end