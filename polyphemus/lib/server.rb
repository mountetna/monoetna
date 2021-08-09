require "json"
require "rack"

require_relative "./polyphemus"
require_relative "./polyphemus/controllers/configuration_controller"
require_relative "./polyphemus/controllers/job_controller"
require_relative "./polyphemus/controllers/ingest_controller"

class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
      application.setup_sequel
    end

    # Return the app host configuration values
    get "/configuration", as: :configuration, action: "configuration#action", auth: { noauth: true }

    # Submit a job for a project
    post "/:project_name/job", action: "job#submit", auth: { user: { is_admin?: :project_name } }

    get "/" do erb_view(:client) end
  end
end
