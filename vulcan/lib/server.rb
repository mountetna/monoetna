require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/browse_controller'
require_relative './server/controllers/workflows_controller'
require_relative './server/controllers/steps_controller'

class Vulcan
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    # root path
    get '/', as: :root do
      erb_view(:client)
    end

    get 'api/workflows', action: 'workflows#fetch', as: :workflows_view, auth: { user: { active?: true } }
    get 'api/workflows/:workflow_name/steps', action: 'steps#fetch', as: :steps_view, auth: { user: { active?: true } }

    with auth: { user: { can_view?: :project_name } } do

      # remaining view routes are parsed by the client and must also be set there
      get ':project_name', as: :project do
        erb_view(:client)
      end

      get ':project_name/*view_path', as: :project_view do
        erb_view(:client)
      end
    end

    def initialize
      super
      # application.setup_db
    end
  end
end
