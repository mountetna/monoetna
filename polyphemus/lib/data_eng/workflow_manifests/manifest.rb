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
          return "Secrets for #{self.as_json[:name]} must be one of: #{self.as_json[:secrets].join(', ')}"
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

      def runtime_params
        as_json[:runtime_params]
      end

      def workflow_path
        as_json[:workflow_path]
      end

      def validate_runtime_config(params)
        errors = []
        params.each do |param, value|
          unless runtime_params.has_key?(param)
            errors.push("no such param #{param}")
            next
          end
          value_opts = runtime_params[param]
          case value_opts[:type]
          when 'select'
            options = value_opts[:values].map { |v| v[:value] }
            unless options.include?(value)
              errors.push("#{param} must be in: #{options.join(', ')}")
            end
          when 'options'
            unless value.is_a?(String) || value_opts.has_key?(:default) && value == value_opts[:default]
              errors.push("#{param} must be a comma-separated string")
            end
          when 'string'
            unless value.is_a?(String) || value_opts.has_key?(:default) && value == value_opts[:default]
              errors.push("#{param} must be a string")
            end
          when 'integer'
            unless value.is_a?(Integer) || value_opts.has_key?(:default) && value == value_opts[:default]
              errors.push("#{param} must be an integer")
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
require_relative './ingestion'
