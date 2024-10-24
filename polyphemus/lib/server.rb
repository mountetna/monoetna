require "json"
require "rack"

require_relative "./polyphemus"
require_relative "./polyphemus/controllers/configuration_controller"
require_relative './polyphemus/controllers/etl_controller'

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

    get '/api/etl/:project_name/configs', action: 'etl#list', auth: { user: { can_edit?: :project_name } }
    post '/api/etl/configs', action: 'etl#list_all', auth: { user: { is_supereditor?: :true } }
    get '/api/etl/jobs', action: 'etl#jobs'
    post '/api/etl/:project_name/update/:config_id', action: 'etl#update', auth: { user: { can_edit?: :project_name } }, log_redact_keys: [ :secrets ]
    post '/api/etl/:project_name/output/:config_id', action: 'etl#add_output', auth: { user: { can_edit?: :project_name } }, log_redact_keys: [ :secrets ]
    post '/api/etl/:project_name/create', action: 'etl#create', auth: { user: { can_edit?: :project_name } }
    get '/api/etl/:project_name/revisions/:config_id', action: 'etl#revisions', auth: { user: { can_edit?: :project_name } }
    get '/api/etl/:project_name/output/:config_id', action: 'etl#output', auth: { user: { can_edit?: :project_name } }

    get '/' do erb_view(:client) end

    get '/:project_name' do erb_view(:client) end
  end
end
