module Etna
  module Codegen
    module Ruby
      class HasSchemasContext
        def schemas_file
          project_path = @key[:project_path]
          project_path.project.schemas_file(project_path.path)
        end

        def resolve_maybe_reference(maybe_ref, cls)

        end
      end

      class EtnaClientFile < File
        include HasFileContext

        def client_class
          sub_context(EtnaClient, identifier: app_name.capitalize)
        end

        def app_name
          @key[:app_name]
        end
      end

      class EtnaClient < Definition
        include HasFileContext
        include Templated
        render_with_erb 'etna_client.rb.erb'

        def servers
          schema.dig('servers')
        end

        def description
          schema.dig('info', 'description') || ''
        end

        def method(method_name)
          sub_context(method_name: method_name) { |key| EtnaClientPathMethod.new(key, flattened_schema_methods[method_name]) }
        end

        def all_methods
          flattened_schema_methods.keys.map { |k| method(k) }
        end

        private

        def flattened_schema_methods
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
        include HasFileContext
        include Templated
        render_with_erb 'etna_client_path_method.rb.erb'

        def initialize(key, method_definition)
          super(key)
          @method_definition = method_definition
        end

        def params_name
          @params_name ||= @key[:method_name] + "_request"
        end

        def params_default
          @params_default ||= imported_identifier(
              @key[:project]
          ) + ".new"
        end
      end
    end
  end
end
