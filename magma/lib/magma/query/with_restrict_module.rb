class Magma
  module WithRestrictModule
    def restrict_constraints_for_model_alias(model, alias_name)
      return {1 => 1} unless @restrict && model_restricted_attribute(model)

      Sequel.negate(
        Sequel.qualify(alias_name, model_restricted_attribute(model).column_name) => true
      )
    end

    def model_restricted_attribute(model)
      model.attributes[:restricted]
    end
  end
end