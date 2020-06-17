require_relative './server/controllers/archimedes_controller'

class Archimedes
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    with auth: { user: { can_view?: :project_name } } do
      # archimedes_controller.rb
      post ':project_name', action: 'archimedes#consignment', as: :consignment
    end

    def initialize(config)
      super
      application.setup_db
      application.setup_magma
    end
  end
end
