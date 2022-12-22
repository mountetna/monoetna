class Magma
  class FloatAttribute < ColumnAttribute
    def database_type
      Float
    end

    def revision_to_loader(record_name, new_value)
      [ name, new_value.to_f ]
    end

    def predicate_class
      Magma::NumberPredicate
    end
  end
end
