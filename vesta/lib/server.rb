require "json"
require "rack"

require_relative "./vesta"
require_relative "./vesta/controllers/project_controller"
require_relative "./vesta/controllers/stats_controller"
require_relative "./vesta/controllers/contact_controller"

class Vesta
  class Server < Etna::Server
    get "/stats", action: "stats#global_stats", auth: { noauth: true }
    get "/stats/projects", action: "stats#project_stats", auth: { noauth: true }

    get "/projects", action: "project#projects", auth: { noauth: true }

    post "/contact", action: "contact#send_contact_email", auth: { noauth: true }

    def initialize
      super
      application.setup_db
      application.setup_sequel
    end
  end
end
