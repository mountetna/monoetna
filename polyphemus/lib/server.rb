class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
    end

      # Return the app host configuration values
      get '/configuration', as: :configuration, action: 'configuration#action', auth: { user: { is_superuser?: :project_name } }

      get '/' do
        [ 200, {}, [ 'Polyphemus is available.' ] ]
      end
  end
end
