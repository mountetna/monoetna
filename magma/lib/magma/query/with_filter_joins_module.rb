class Magma
  module WithFilterJoinsModule
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