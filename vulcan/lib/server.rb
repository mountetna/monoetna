require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/browse_controller'
require_relative './server/controllers/workflows_controller'
require_relative './server/controllers/data_controller'
require_relative './server/controllers/sessions_controller'
require_relative './server/models/session'
require_relative './server/models/workflow'

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

    get 'api/:project_name/workflows', action: 'workflows#fetch', as: :workflows_view, auth: { user: { can_view?: :project_name } }
    get 'api/:project_name/data/:cell_hash/:data_filename', action: 'data#fetch', as: :data_view, match_ext: true, auth: { user: { can_view?: :project_name } }
    post 'api/:project_name/session/:workflow_name', action: 'sessions#submit', as: :submit_view, match_ext: true, auth: { user: { can_view?: :project_name } }

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
