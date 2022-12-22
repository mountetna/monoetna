class Magma
  class BooleanAttribute < ColumnAttribute
    def database_type
      TrueClass
    end

    def predicate_class
      Magma::BooleanPredicate
    end
  end
end
