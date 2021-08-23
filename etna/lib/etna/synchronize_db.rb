module Etna
  class SynchronizeDb
    def initialize(app)
      @app = app
    end

    def call(env)
      # Do a coarse checkout of the connection
      Etna::Application.instance.db.synchronize do
        @app.call(env)
      end
    end
  end
end
