require_relative 'templates'
require 'json'



module Etna
  module Codegen
    module Ruby
      class RequiresSet < Context
        include Templated
        render_with_erb 'requires.rb.erb'

        def add_require(project_path:, **key)
          path = @key[:project_path].relative_to(project_path)
          sub_context(path) { path }
        end
      end

      module HasFileContext
        def imported_identifier(scopes:, project_path:, identifier:, **key)
          @key[:file].requires.add_require(project_path: project_path)
          "#{relative_scope(scopes: scopes)}::#{identifier}"
        end
      end

      class File < Context
        include Templated
        render_with_erb 'ruby_file.rb.erb'

        def initialize(key)
          super({file: self}.update(key))
        end

        def requires
          @requires ||= RequiresSet.new(@key)
        end

        def relative_scope(scopes:, **key)
          self_scopes = @key[:scopes].values.flatten
          other_scopes = scopes.values.flatten

          common_scopes = 0
          self_scopes.zip(other_scopes).each do |self_scope, other_scope|
            if self_scope != other_scope
              break
            end

            common_scopes += 1
          end

          other_scopes.slice(common_scopes..-1).join('::')
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
      end

      class Definition < Context
        include HasFileContext

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
        include HasFileContext

        def method_name
          @key[:method_name]
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
