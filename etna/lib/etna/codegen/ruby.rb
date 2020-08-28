require_relative 'templates'

module Etna
  module Codegen
    module Ruby
      class File < Context
        include Templated
        render_with_erb 'ruby_file.erb'

        def relative_requires
          @relative_requires ||= Set.new
        end

        def requires
          @requires ||= Set.new
        end

        def add_require(project_path, *exported)
          exported
        end

        def module_scopes
          @module_scopes ||= begin
            scopes = @key[:scopes] || []
            result = []
            parent = ""
            scopes.each do |key, scope|
              scope.split('::').each do |m|
                parent = "#{parent}::#{m}"
                result << "#{key} #{parent}"
              end
            end

            result
          end
        end

        protected

        def identifiers
          @values.values
        end
      end

      class Definition < Context
        def identifier
          @key[:identifier]
        end

        def commented(text)
          text.split('\n').map do |t|
            "# #{t}"
          end.join('\n')
        end
      end

      class MethodDefinition < Definition
        def identifier
          @key[:method_name]
        end
      end

      class EtnaClientFile < File
        def client_class
          sub_context(EtnaClient, identifier: app_name.capitalize)
        end

        def app_name
          @key[:app_name]
        end
      end

      class EtnaClient < Definition
        include Templated
        render_with_erb 'etna_client.rb.erb'

        def servers
          schema.dig('servers')
        end

        def description
          schema.dig('info', 'description') || ''
        end

        def method(method_name)
          sub_context(method_name: method_name) { |key| EtnaClientPathMethod.new(key, flattened_methods[method_name]) }
        end

        def all_methods
          flattened_methods.keys.map { |k| method(k) }
        end

        private

        def methods
          @values.values
        end

        def flattened_methods
          result = {}
          @flattened_methods ||= schema.dig('paths').each do |path, http_methods|
            http_methods.each do |http_method, v|
              method_name = v['x-method-name']
              if method_name.nil?
                raise "Api #{@key[:app_name]}[#{http_method} #{path}] is lacking x-method-name, cannot generate client"
              end

              if result.include?(method_name)
                raise "Api #{@key[:app_name]}[#{http_method} #{path}] has duplicate definitions for x-method-name #{method_name}"
              end

              result[method_name] = v
            end
          end

          result
        end

        def schema
          @key[:swagger_file].contents
        end
      end

      class EtnaClientPathMethod < MethodDefinition
        include Templated
        render_with_erb 'etna_client_path_method.rb.erb'

        def initialize(key, method_definition)
          super(key)
          @method_definition = method_definition
        end
      end
    end
  end
end
