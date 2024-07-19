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

    # Vulcan V2 endpoints

    # CRUD Repo
    post 'api/v2/repo/create', action: 'vulcan_v2#create_repo', auth: { user: { is_admin?: :project_name }}
    get 'api/v2/:project_name/repo/', action: 'vulcan_v2#list_repos', auth: { user: { is_admin?: :project_name }}

    # CRUD Workflow
    post 'api/v2/workflow/publish', action: 'vulcan_v2#publish_workflow'
    get 'api/v2/:project_name/workflows/', action: 'vulcan_v2#list_workflows'

    # CRUD Workspace
    post 'api/v2/:project_name/workspace/create', action: 'vulcan_v2#create_workspace'
    get 'api/v2/:project_name/workspace', action: 'vulcan_v2#list_workspaces'
    get 'api/v2/:project_name/workspace/:workspace_id', action: 'vulcan_v2#get_workspace'

    # Run API
    post 'api/v2/:project_name/workspace/:workspace_id/run', action: 'vulcan_v2#run_workflow'
    get'api/v2/:project_name/workspace/:workspace_id/run/:run_id', action: 'vulcan_v2#get_workflow_status'

    # File API
    post 'api/v2/:project_name/workspace/:workspace_id/file/write', action: 'vulcan_v2#write_files'
    post 'api/v2/:project_name/workspace/:workspace_id/file/read', action: 'vulcan_v2#read_files'


    # Vulcan V1 endpoints - to remove
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



    end

    # root path
    get '/', as: :root do erb_view(:client) end

    def initialize
      super
      application.setup_db
      application.setup_ssh
    end
  end
end
