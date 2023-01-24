class Magma
  module WithDateTimePredicateMethods
    def extract table, identity, is_all
      Magma::Answer.new(table.first[column_name]&.iso8601)
    end
  end
end