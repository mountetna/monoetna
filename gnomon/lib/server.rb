require "json"
require "rack"

require_relative "./gnomon"
require_relative "./gnomon/controllers/gnomon_controller"

class Gnomon
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
      application.setup_sequel
    end

    get '/' do erb_view(:client) end

    get '/:project_name' do erb_view(:client) end

    get '/:project_name/*etc' do erb_view(:client) end
  end
end
