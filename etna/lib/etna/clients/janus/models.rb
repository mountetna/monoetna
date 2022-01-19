require 'ostruct'
require_relative '../../json_serializable_struct'
require_relative '../base_client'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Janus < Etna::Clients::BaseClient
      class GetProjectRequest < Struct.new(:project_name, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end
      end

      class AddProjectRequest < Struct.new(:project_name, :project_name_full, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end
      end

      class AddUserRequest < Struct.new(:project_name, :email, :name, :role, :privileged, :affiliation, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({privileged: false, affiliation: ''}.update(params))
        end
      end

      class UpdatePermissionRequest < Struct.new(:project_name, :email, :role, :privileged, :affiliation, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({privileged: false, affiliation: ''}.update(params))
        end
      end

      class RefreshTokenRequest
        def map
          []
        end
      end

      class ViewerTokenRequest
        def map
          []
        end
      end

      class ValidateTaskTokenRequest
        def map
          []
        end
      end

      class HtmlResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def html
          @raw
        end
      end

      class GetProjectResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def project
          Project.new(@raw[:project])
        end
      end

      class GetProjectsResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def projects
          @raw[:projects].map do |project|
            Project.new(project)
          end
        end
      end

      class Project
        def initialize(raw = '')
          @raw = raw
        end

        def project_name
          @raw[:project_name]
        end

        def project_name_full
          @raw[:project_name_full]
        end

        def permissions
          @raw[:permissions] || []
        end

        def resource
          !!@raw[:resource]
        end
      end

      class TokenResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def token
          @raw
        end
      end
    end
  end
end
