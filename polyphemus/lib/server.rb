require 'json'
require 'rack'

require_relative './polyphemus'
require_relative './polyphemus/controllers/configuration_controller'

class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
    end

    # Return the app host configuration values
    get '/configuration', as: :configuration, action: 'configuration#action', auth: { noauth: true }

    get '/' do
      [ 200, {}, [ 'Polyphemus is available.' ] ]
    end
  end
end
