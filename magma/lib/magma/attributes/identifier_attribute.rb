class Magma
  class IdentifierAttribute < StringAttribute
    def initialize(opts = {})
      @primary_key = opts.delete(:primary_key)
      super(opts.merge(unique: true))
    end

    private

    def after_magma_model_set
      @magma_model.identity = self
      @magma_model.order(column_name) unless @magma_model.order
    end

    class Validation < Magma::Validation::Attribute::BaseAttributeValidation
      def validate(value, &block)
        gnomon_mode = @flags[Magma::Flags::GnomonMode::NAME]
        return if gnomon_mode.nil? or gnomon_mode == Magma::Flags::GnomonMode::NONE
        if gnomon_mode == Magma::Flags::GnomonMode::IDENTIFIER
          yield(
            "The identifier '#{value}' has not been assigned in gnomon."
          ) unless Magma.instance.db[:identifiers].where(identifier: value).any?
        end
      end
    end
  end
end
