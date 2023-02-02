class Magma
  module WithSubqueriesModule
    def apply_subqueries(query)
      return query if @subqueries.nil? || @subqueries.empty?

      @subqueries.each do |subquery|
        query = subquery.apply(query)
      end

      query
    end
  end
end