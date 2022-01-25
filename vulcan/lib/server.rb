require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/browse_controller'
require_relative './server/controllers/workflows_controller'
require_relative './server/controllers/data_controller'
require_relative './server/controllers/sessions_controller'

class Vulcan
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    with auth: { user: { active?: true, has_flag?: 'vulcan' } } do
      get 'api/workflows', action: 'workflows#fetch', as: :workflows_view
      get 'api/:project_name/workflows', action: 'workflows#fetch_for_project'
      get 'api/:project_name/data/:cell_hash/:data_filename', action: 'data#fetch', as: :data_view, match_ext: true
      post 'api/:project_name/session/:workflow_name/status', action: 'sessions#status', as: :status_view, match_ext: true
      post 'api/:project_name/session/:workflow_name', action: 'sessions#submit', as: :submit_view, match_ext: true
      get 'api/:project_name/workflows/:workflow_name/session/:key', action: 'sessions#fetch'
      post 'api/:project_name/workflows/:workflow_name/from_query', action: 'sessions#create_from_query'

      # remaining view routes are parsed by the client and must also be set there
      get '/:project_name/workflow/:session_key', action: 'sessions#render'

      get '/:project_name/workflow', as: :workflow do
        erb_view(:client)
      end

      get '/:project_name/workflow/*view_path', as: :workflow_view do
        erb_view(:client)
      end
    end

    # root path
    get '/', as: :root do erb_view(:client) end
    get '/:project_name', as: :root do erb_view(:client) end

    def initialize
      super
      application.setup_db
    end
  end
end
