class Magma
  class MatchPredicate < Magma::ColumnPredicate
    verb '::type' do
      child String

      extract do |table, identity|
        Magma::Answer.new(table.first[column_name]['type'])
      end
    end

    verb '::value' do
      child String

      extract do |table, identity|
        Magma::Answer.new(table.first[column_name]['value'])
      end
    end

    verb nil do
      child String
    end

    def select(incoming_alias_name=nil, incoming_attribute=nil)
      [ Sequel[alias_name][@column_name].as(column_name) ]
    end
  end
end
