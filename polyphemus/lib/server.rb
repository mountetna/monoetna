require 'json'
require 'rack'

require_relative './polyphemus'
require_relative './polyphemus/controllers/configuration_controller'
require_relative './polyphemus/controllers/job_controller'
require_relative './polyphemus/controllers/ingest_controller'

class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
      application.setup_sequel
    end

    # Return the app host configuration values
    get '/configuration', as: :configuration, action: 'configuration#action', auth: { noauth: true }

    # Submit a job for a project
    post '/:project_name/job', action: 'job#submit', auth: { user: { is_admin?: :project_name } }

    # Get files associated with a directory, for ingest
    get '/:project_name/ingest/list/:ingest_host/*folder_path', action: 'ingest#list_dir', auth: { user: { is_admin?: :project_name } }, match_ext: true
    # Update the status of ingest files associated with the directory
    post '/:project_name/ingest/update/:ingest_host/*folder_path', action: 'ingest#update', auth: { user: { is_admin?: :project_name } }, match_ext: true

    get '/' do erb_view(:client) end
  end
end
