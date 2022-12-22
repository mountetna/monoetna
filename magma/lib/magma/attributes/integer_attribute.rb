class Magma
  class IntegerAttribute < ColumnAttribute
    def database_type
      Integer
    end

    def revision_to_loader(record_name, new_value)
      [ name, new_value.to_i ]
    end

    def predicate_class
      Magma::NumberPredicate
    end
  end
end
