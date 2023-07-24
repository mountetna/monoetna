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
        if flag_value.nil? or flag_value == gnomon_mode[:none]
          # Some identifier attributes have additional validations, so we call super() to run those validations first
          super if @attribute.validation
          return
        end
        begin
        # Identifier mode - we just care about the existence of the identifier in the identifier table
          if flag_value == gnomon_mode[:identifier]
            identifier = Magma::Gnomon::Identifier.where(project_name: @model.project_name.to_s, identifier: value).first
            raise "The identifier '#{value}' has not been assigned in Gnomon." unless identifier

          # Pattern mode - We must make sure the identifier conforms to a grammar
          elsif flag_value == gnomon_mode[:pattern]
            grammar = Magma::Gnomon::Grammar.for_project(@model.project_name.to_s)
            raise "Grammar not found, identifier is: '#{value}'" unless grammar

            rule = grammar.parser.rules[@model.model_name.to_s]
            raise "No such rule #{@model.model_name} for #{@model.project_name}" unless rule
            raise "The identifier '#{value}' does not conform to a grammar in Gnomon." unless rule.valid?(value)
          end
        rescue Exception => e
          yield "#{e.message}"
        end
      end

    end
  end
end
