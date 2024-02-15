require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/browse_controller'
require_relative './server/controllers/workflows_controller'
require_relative './server/controllers/figure_controller'
require_relative './server/controllers/data_controller'
require_relative './server/controllers/sessions_controller'
require_relative './server/controllers/vulcan_v2_controller'

class Vulcan
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    get 'api/workflows', action: 'workflows#fetch', as: :workflows_view, auth: { user: { active?: true } }

    with auth: { user: { can_view?: :project_name } } do
      get 'api/:project_name/data/:cell_hash/:data_filename', action: 'data#fetch', as: :data_view, match_ext: true
      post 'api/:project_name/session/:workflow_name/status', action: 'sessions#status', as: :status_view, match_ext: true
      post 'api/:project_name/session/:workflow_name', action: 'sessions#submit', as: :submit_view, match_ext: true
      get 'api/:project_name/workflows', action: 'workflows#fetch_for_project'

      get 'api/:project_name/figures', action: 'figure#fetch'
      get 'api/:project_name/figure/:figure_id', action: 'figure#get'
      get 'api/:project_name/figure/:figure_id/revisions', action: 'figure#revisions'
      post 'api/:project_name/figure/create', action: 'figure#create'
      post 'api/:project_name/figure/:figure_id/update', action: 'figure#update'
      delete 'api/:project_name/figure/:figure_id', action: 'figure#delete'

      # remaining view routes are parsed by the client and must be set there
      get '/:project_name', as: :project_root do erb_view(:client) end
      get '/:project_name/*client_path', as: :client_view do erb_view(:client) end

      # vulcan v2
      post 'api/:project_name/workflow/list', action: 'vulcan_v2#list'
      post 'api/:project_name/workflow/init', action: 'vulcan_v2#init'

    end

    # root path
    get '/', as: :root do erb_view(:client) end

    def initialize
      super
      application.setup_db
    end
  end
end
