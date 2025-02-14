require "json"
require "rack"

require_relative "./polyphemus"
require_relative "./polyphemus/controllers/configuration_controller"
require_relative './polyphemus/controllers/workflow_controller'

class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
      application.setup_sequel
      application.setup_ssh
    end

    # Return the app host configuration values
    get "/configuration", as: :configuration, action: "configuration#action", auth: { noauth: true }

    get "/:project_name/ingest/hosts", action: "ingest#list_hosts", auth: { user: { is_admin?: :project_name } }
    # Get files associated with a directory, for ingest
    get "/:project_name/ingest/list/:ingest_host/*folder_path", action: "ingest#list_dir", auth: { user: { is_admin?: :project_name } }, match_ext: true
    # Add files in the folder to the ingest queue
    post "/:project_name/ingest/enqueue/:ingest_host/*folder_path", action: "ingest#enqueue", auth: { user: { is_admin?: :project_name } }, match_ext: true

    get '/' do erb_view(:client) end
    get '/:project_name' do erb_view(:client) end

    # Polyphemus API V2

    post '/api/workflows/:project_name/create', action: 'workflow#create', auth: { user: { can_edit?: :project_name } }
    post '/api/workflows/:project_name/update/:config_id', action: 'workflow#update', auth: { user: { can_edit?: :project_name } }, log_redact_keys: [ :secrets ]

    get '/api/workflows/:project_name/configs/:config_id', action: 'workflow#list_by_id', auth: { user: { can_edit?: :project_name } }
    post '/api/workflows/configs', action: 'workflow#list_all', auth: { user: { is_supereditor?: :true } }

    get '/api/workflows', action: 'workflow#get_workflows'
    get '/api/workflows/:project_name/revisions/:config_id', action: 'workflow#revisions', auth: { user: { can_view?: :project_name } }

    get '/api/workflows/:project_name/runs/:config_id', action: 'workflow#runs', auth: { user: { can_edit?: :project_name } }
    get '/api/workflows/:project_name/run/:run_id', action: 'workflow#get_run', auth: { user: { can_view?: :project_name } }
    get '/api/workflows/:project_name/run/output/:run_id', action: 'workflow#run_output', auth: { user: { can_view?: :project_name } }
    post '/api/workflows/:project_name/run/update/:run_id', action: 'workflow#update_run', auth: { user: { can_edit?: :project_name } }
    post '/api/workflows/:project_name/run/previous/:config_id', action: 'workflow#get_previous_state', auth: { user: { can_edit?: :project_name } }

    get '/api/workflows/:project_name/runtime_configs/:config_id', action: 'workflow#get_runtime_config', auth: { user: { can_view?: :project_name } }
    post '/api/workflows/:project_name/runtime_configs/update/:config_id', action: 'workflow#update_runtime_config', auth: { user: { can_edit?: :project_name } }

    post '/api/workflows/:project_name/runtime_configs/run_once/:config_id', action: 'workflow#run_once', auth: { user: { can_edit?: :project_name } }
    get '/api/workflows/:project_name/status', action: 'workflow#status', auth: { user: { can_edit?: :project_name } }
  end
end
