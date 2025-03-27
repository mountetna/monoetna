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

    # Workflow API
    post 'api/v2/:project_name/workflows/create', action: 'vulcan_v2#create_workflow', auth: { user: { is_admin?: :project_name }}
    get 'api/v2/:project_name/workflows/', action: 'vulcan_v2#list_workflows', auth: { user: { can_edit?: :project_name }}

    # Workspace API
    post 'api/v2/:project_name/workspace/create', action: 'vulcan_v2#create_workspace', auth: { user: { can_edit?: :project_name }}
    get 'api/v2/:project_name/workspace', action: 'vulcan_v2#list_workspaces',  auth: { user: { can_edit?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id', action: 'vulcan_v2#get_workspace',  auth: { user: { can_edit?: :project_name }}

    # Config API
    post 'api/v2/:project_name/workspace/:workspace_id/config', action: 'vulcan_v2#save_config', auth: { user: { can_edit?: :project_name }}
    get'api/v2/:project_name/workspace/:workspace_id/config', action: 'vulcan_v2#get_config', auth: { user: { can_edit?: :project_name }}

    # Run API
    post 'api/v2/:project_name/workspace/:workspace_id/run/:config_id', action: 'vulcan_v2#run_workflow', auth: { user: { can_edit?: :project_name }}
    get'api/v2/:project_name/workspace/:workspace_id/run/:run_id', action: 'vulcan_v2#get_workflow_status',  auth: { user: { can_edit?: :project_name }}

    # DAG endpoint
    get'api/v2/:project_name/workspace/:workspace_id/dag', action: 'vulcan_v2#get_dag',  auth: { user: { can_edit?: :project_name }}

    # File API
    post 'api/v2/:project_name/workspace/:workspace_id/file/write', action: 'vulcan_v2#write_files', auth: { user: { can_edit?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/file/read', action: 'vulcan_v2#read_files', auth: { user: { can_edit?: :project_name }}
    get'api/v2/:project_name/workspace/:workspace_id/file/', action: 'vulcan_v2#get_files', auth: { user: { can_edit?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/image/read', action: 'vulcan_v2#read_image', auth: { user: { can_edit?: :project_name }}

    # Is running endpoint
    get'api/v2/:project_name/workspace/:workspace_id/running', action: 'vulcan_v2#is_running',  auth: { user: { can_edit?: :project_name }}
  
    # Vulcan V1 endpoints - to eventually remove
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
