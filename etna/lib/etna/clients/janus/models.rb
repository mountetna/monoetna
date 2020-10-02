require 'ostruct'
require_relative '../../json_serializable_struct'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Janus
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

      class WhoAmIRequest
        def map
          []
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

      class HtmlResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def html
          @raw
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

      class UserResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def user
          User.new(@raw)
        end
      end

      class User
        def initialize(raw = '')
          @raw = raw
        end

        def email
          @raw['email']
        end

        def first
          @raw['first']
        end

        def last
          @raw['last']
        end
      end
    end
  end
end