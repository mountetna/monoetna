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

        gnomon_mode = Magma::Flags::GNOMON_MODE
        flag_value = @flags[gnomon_mode[:name]]

        # We skip validation if we have the gnomon flag or if it's value is set to none
        return if flag_value.nil? or flag_value == gnomon_mode[:none]

        # Identifier mode - we just care about the existence of the identifier in the identifier table
        if flag_value == gnomon_mode[:identifier]
          return if Magma.instance.db[:identifiers].where(project_name: @model.project_name.to_s, identifier: value).any?
          yield "The identifier '#{value}' has not been assigned in Gnomon."

          # Pattern mode - We must make sure the identifier conforms to a grammar
        elsif flag_value == gnomon_mode[:pattern]
          identifier = Magma.instance.db[:identifiers].where(project_name: @model.project_name.to_s, identifier: value).first
          grammar =  Magma.instance.db[:grammars].where(id: identifier[:grammar_id]).first
          gnomon_validation = Magma::Gnomon::Validation.new( Magma::Gnomon::Grammar::Parser.new(grammar[:config]) )

          return if gnomon_validation.valid?
          yield "The identifier '#{value}' does not conform to a grammar in Gnomon."
        end
      end
    end
  end
end
