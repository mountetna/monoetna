class Magma
  class Validation
    class Project < Magma::Validation::Model
      def initialize(model, validator)
        super
      end

      def project_name
        @project_name ||= @model.select_map( @model.identity.column_name.to_sym ).first
      end

      def grammar_project_name
        @grammar_project_name ||= @validator.grammar&.token_project_name
      end

      def self.skip?(model)
        model.model_name != :project
      end

      def valid_names
        return @valid_names if @valid_names

        @valid_names = [
          project_name,
          @model.project_name,
          @model.project_name.upcase,
          @validator.grammar&.config&.dig("tokens","PROJ","values")&.keys&.first,
          grammar_project_name
        ].compact.map(&:to_s)
      end

      def validate(record_name, document)
        return if !project_name && !grammar_project_name

        unless valid_names.include?(record_name)
          yield "Project name must match one of: #{valid_names.join(', ')}"
        end
      end
    end
  end
end
