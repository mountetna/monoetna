class Polyphemus
  class Server < Etna::Server
    def initialize
      super
      # application.load_models
      application.setup_db
    end
  end
end
