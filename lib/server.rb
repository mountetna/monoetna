# routes.rb
# This file initiates Metis and sets the routing of the http requests.

class Metis
  class Server < Etna::Server
    get '/', 'client#index'
    get '/user', 'client#user'

    post '/check', 'user_log#check_log'

    post '/upload/authorize', 'upload#authorize'
    post '/upload/start', 'upload#start'
    post '/upload/blob', 'upload#blob'
    post '/upload/pause', 'upload#pause'
    post '/upload/cancel', 'upload#cancel'
    post '/upload/remove-file', 'upload#remove_file'
    post '/upload/remove-failed', 'upload#remove_failed'
    post '/upload/recover', 'upload#recover'

    post '/file', 'download#file'

    def initialize(config)
      super
      application.load_models
    end
  end
end
