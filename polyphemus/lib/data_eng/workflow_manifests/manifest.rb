class Polyphemus
  class WorkflowManifest
    class << self
      # Allows us to dynamically add subclasses to the list
      # and create methods on subclasses
      def inherited(subclass)
        @list ||= []
        @list << subclass
      end
      attr_reader :list

      def from_workflow_name(workflow_name)
        list.find do |workflow|
          workflow.as_json[:name] == workflow_name
        end
      end

      def validate_secrets(secrets)
        unless (secrets.keys - self.as_json[:secrets]).empty?
          return "Secrets for #{self.as_json[:name]} jobs must be one of: #{self.as_json[:secrets].join(', ')}"
        end
      end

      def validate_config(config_to_validate)
        schema = JSONSchemer.schema(
          JSON.parse(self.as_json[:schema].to_json)
        )
        schema.validate(JSON.parse(config_to_validate.to_json)).to_a
      end
    
      def valid_config?(config)
        schema = JSONSchemer.schema(
        JSON.parse(self.as_json[:schema].to_json)
      )
      schema.valid?(JSON.parse(config.to_json))
      end

      def job_params
        as_json[:params]
      end

      def validate_params(params)
        errors = []
        params.each do |param, value|
          unless job_params.has_key?(param)
            errors.push("no such param #{param}")
            next
          end
          value_opts = job_params[param]
          case value_opts
          when Array
            options = value_opts.map { |v| v[:value] }
            unless options.include?(value)
              errors.push("#{param} must be in: #{options.join(', ')}")
            end
          when Hash
            unless value.is_a?(String)
              errors.push("#{param} must be a comma-separated string")
            end
          when 'string'
            unless value.is_a?(String)
              errors.push("#{param} must be a string")
            end
          when 'boolean'
            unless [true,false].include?(value)
              errors.push("#{param} must be a boolean")
            end
          end
        end
        errors
      end

    end

    attr_reader :request_params, :errors, :token
    def initialize(request_params:, token:)
      @request_params = request_params
      @token = token
      @errors = []
    end

    def valid?
      errors.length == 0
    end

    def validate
      raise "This should be implemented in a subclass."
    end

    def as_json
      raise "This should be implemented in a subclass."
    end

    private

    def require_params(*params)
      missing_params = params.reject{|p| @request_params.key?(p) }
      raise JobError, "request_params missing required param(s): #{missing_params.join(', ')}" unless missing_params.empty?
    end
    alias_method :require_param, :require_params
  end

end

require_relative './metis_linker'
require_relative './redcap_loader'
#require_relative './cat_ingestion'
