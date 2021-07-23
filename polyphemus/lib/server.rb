require "json"
require "rack"

require_relative "./polyphemus"
require_relative "./polyphemus/controllers/configuration_controller"
require_relative "./polyphemus/controllers/job_controller"
require_relative "./polyphemus/controllers/ingest_controller"
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

    # Submit a job for a project
    post "/:project_name/job", action: "job#submit", auth: { user: { is_admin?: :project_name } }

    get "/:project_name/ingest/hosts", action: "ingest#list_hosts", auth: { user: { is_admin?: :project_name } }
    # Get files associated with a directory, for ingest
    get "/:project_name/ingest/list/:ingest_host/*folder_path", action: "ingest#list_dir", auth: { user: { is_admin?: :project_name } }, match_ext: true
    # Add files in the folder to the ingest queue
    post "/:project_name/ingest/enqueue/:ingest_host/*folder_path", action: "ingest#enqueue", auth: { user: { is_admin?: :project_name } }, match_ext: true

    get '/api/:project_name/etl/configs', action: 'etl#list_configs', auth: { user: { can_edit?: :project_name } }
    get '/api/:project_name/etl/jobs', action: 'etl#list_jobs', auth: { user: { can_edit?: :project_name } }

    get '/' do erb_view(:client) end

    get '/:project_name' do erb_view(:client) end
  end
end
