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

    get '/api/etl/:project_name/configs', action: 'workflow#list', auth: { user: { can_edit?: :project_name } }
    #post '/api/etl/configs', action: 'etl#list_all', auth: { user: { is_supereditor?: :true } }
    get '/api/etl/jobs', action: 'etl#jobs'
    # Change format of this
    post '/api/etl/:project_name/output/:config_id', action: 'workflow#add_output', auth: { user: { can_edit?: :project_name } }, log_redact_keys: [ :secrets ]
    post '/api/etl/:project_name/create', action: 'workflow#create', auth: { user: { can_edit?: :project_name } }
    get '/api/etl/:project_name/revisions/:config_id', action: 'workflow#revisions', auth: { user: { can_edit?: :project_name } }
    get '/api/etl/:project_name/output/:config_id', action: 'workflow#output', auth: { user: { can_edit?: :project_name } }

    get '/' do erb_view(:client) end

    get '/:project_name' do erb_view(:client) end

    # Polyphemus API V2
    post '/api/etl/:project_name/create', action: 'workflow#create', auth: { user: { can_edit?: :project_name } }
    post '/api/etl/:project_name/update/:config_id', action: 'workflow#update', auth: { user: { can_edit?: :project_name } }, log_redact_keys: [ :secrets ]


    get "/api/elt/:project_name/configs/:config_id", action: "workflow#get_config", auth: { user: { can_edit?: :project_name } }
    get "/api/etl/:project_name/:config_id/run/:run_id", action: "workflow#get_workflow_state", auth: { user: { can_edit?: :project_name } }
    post "/api/etl/:project_name/:config_id/run/:run_id", action: "workflow#update_workflow_state", auth: { user: { can_edit?: :project_name } }
    post "/api/etl/:project_name/:config_id/metadata", action: "workflow#write_run_metadata", auth: { user: { can_edit?: :project_name } }
  end
end
