require_relative '../templates'
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

    end
  end
end
