class Magma
  class StringAttribute < ColumnAttribute
    def database_type
      String
    end

    def predicate_class
      Magma::StringPredicate
    end
  end
end
