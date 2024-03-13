class Magma
  class StringAttribute < Attribute
    def database_type
      String
    end

    def revision_to_loader(record_name, new_value)
      [ name, new_value&.to_s ]
    end
  end
end
