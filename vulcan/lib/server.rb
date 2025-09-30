require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/vulcan_v2_controller'

class Vulcan
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    # Workflow API
    post 'api/v2/:project_name/workflows/create', action: 'vulcan_v2#create_workflow', auth: { user: { is_superuser?: :project_name }}
    get 'api/v2/:project_name/workflows/', action: 'vulcan_v2#list_workflows', auth: { user: { can_view?: :project_name }}

    # Workspace API
    post 'api/v2/:project_name/workspace/create', action: 'vulcan_v2#create_workspace', auth: { user: { can_view?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/update', action: 'vulcan_v2#update_workspace', auth: { user: { can_view?: :project_name }}
    delete 'api/v2/:project_name/workspace/:workspace_id', action: 'vulcan_v2#delete_workspace', auth: { user: { is_admin?: :project_name }}
    get 'api/v2/:project_name/workspace', action: 'vulcan_v2#list_workspaces',  auth: { user: { can_view?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id', action: 'vulcan_v2#get_workspace',  auth: { user: { can_view?: :project_name }}

    # Config API
    post 'api/v2/:project_name/workspace/:workspace_id/config', action: 'vulcan_v2#save_config', auth: { user: { can_view?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id/config/:config_id', action: 'vulcan_v2#get_config', auth: { user: { can_view?: :project_name }}

    # Run API
    post 'api/v2/:project_name/workspace/:workspace_id/run/:config_id', action: 'vulcan_v2#run_workflow', auth: { user: { can_view?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id/run/:run_id', action: 'vulcan_v2#get_workflow_status',  auth: { user: { can_view?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/:run_id/cancel', action: 'vulcan_v2#cancel_workflow', auth: { user: { can_view?: :project_name }}

    # DAG endpoint
    get'api/v2/:project_name/workspace/:workspace_id/dag', action: 'vulcan_v2#get_dag',  auth: { user: { can_view?: :project_name }}

    # File API
    post 'api/v2/:project_name/workspace/:workspace_id/file/write', action: 'vulcan_v2#write_files', auth: { user: { can_view?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/file/read', action: 'vulcan_v2#read_files', auth: { user: { can_view?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id/file/', action: 'vulcan_v2#get_files', auth: { user: { can_view?: :project_name }}
    post 'api/v2/:project_name/workspace/:workspace_id/image/read', action: 'vulcan_v2#read_image', auth: { user: { can_view?: :project_name }}
    get 'api/v2/:project_name/workspace/:workspace_id/file/download/*file_name', action: 'vulcan_v2#download_file', auth: { user: { can_view?: :project_name }}

    # State API
    post 'api/v2/:project_name/workspace/:workspace_id/state', action: 'vulcan_v2#get_state', auth: { user: { can_view?: :project_name }}

    # Is running endpoint
    get 'api/v2/:project_name/workspace/:workspace_id/running', action: 'vulcan_v2#is_running',  auth: { user: { can_view?: :project_name }}
  
    # Cluster latency endpoint
    get 'api/v2/cluster-latency', action: 'vulcan_v2#cluster_latency'
  
    # root path
    get '/', as: :root do erb_view(:client) end

    # Project-specific routes - needed for frontend to work with project URLs
    with auth: { user: { can_view?: :project_name } } do
      get '/:project_name', as: :project_root do erb_view(:client) end
      get '/:project_name/*client_path', as: :client_view do erb_view(:client) end
    end

    def initialize
      super
      application.setup_db
      application.setup_ssh
      application.vulcan_checks
    end
  end
end
