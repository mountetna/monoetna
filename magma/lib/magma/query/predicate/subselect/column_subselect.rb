require_relative '../column'

class Magma
  class ColumnSubselectPredicate < Magma::ColumnPredicate
    def extract table, identity, is_all
      if @verb && @verb.gives?(:extract)
        @verb.do(:extract, table, identity)
      elsif table.first[column_name].nil?
        Magma::NilAnswer.new
      elsif is_all
        nested_reduce_and_apply(
          table.first[column_name],
          0,
          &do_nothing
        )
      else
        Magma::Answer.new(table.first[column_name])
      end
    end

    def select
      [ Magma::SubselectColumn.new(column_name, attribute_column_name) ]
    end

    private

    def do_nothing
      lambda { |d|
        d
      }
    end
  end
end
