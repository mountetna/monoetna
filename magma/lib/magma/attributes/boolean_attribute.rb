class Magma
  class BooleanAttribute < Attribute
    def database_type
      TrueClass
    end

    def revision_to_loader(record_name, new_value)
      v = case new_value
          when 0, "0"
            false
          when 1, "1"
            true
          else
            new_value
          end
      [ name, v ]
    end
  end
end
