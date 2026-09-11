class Magma
  module Gnomon
    class Validation
      attr_reader :errors

      def initialize(grammar, config, comment, project_record_name)
        @grammar = grammar
        @config = config
        @comment = comment
        @project_record_name = project_record_name
        @errors = []
      end

      def valid?
        validations.each do |validation|
          send(validation)
        end

        errors.empty?
      end

      private

      def validations
        [
          :validate_schema,
          :validate_rules,
          :validate_tokens,
          :validate_synonyms,
          :confirm_project_token_change
        ]
      end

      def validate_schema
        schema = JSONSchemer.schema(
          JSON.parse(Magma::Gnomon::Grammar.to_schema.to_json)
        )

        schema_errors = schema.validate(JSON.parse(@grammar.config.to_json))

        @errors += schema_errors.map do |error|
          JSONSchemer::Errors.pretty(error)
        end
      end

      def validate_rules
        @errors += @grammar.rules.errors unless @grammar.rules.valid?
      end

      def validate_tokens
        @errors += @grammar.tokens.errors unless @grammar.tokens.valid?
      end

      def validate_synonyms
        @errors += @grammar.synonyms.errors unless @grammar.synonyms.valid?
      end

      def confirm_project_token_change
        project_token_name = @config.dig('tokens','PROJECT','values')&.first&.first
        if project_token_name
          project_new_hash = Digest::MD5.hexdigest(project_token_name)
          if @project_record_name != project_token_name && !@comment.include?(project_new_hash)
            @errors += ["Add \"#{project_new_hash}\" in comment to confirm intent to change project name to provided PROJECT token value, \"#{project_token_name}\"."]
          end
        end
      end
    end
  end
end
