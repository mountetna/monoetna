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
        gnomon_mode = @flags[Magma::Flags::GNOMON_MODE[:name]]
        return if gnomon_mode.nil? or gnomon_mode == Magma::Flags::GNOMON_MODE[:none]
        if gnomon_mode == Magma::Flags::GNOMON_MODE[:identifier]
          yield(
            "The identifier '#{value}' has not been assigned in Gnomon."
          ) unless Magma.instance.db[:identifiers].where(identifier: value).any?
        end
      end
    end
  end
end
