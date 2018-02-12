require 'json'

require_relative './server/controllers/metis_controller'
require_relative './server/controllers/main_controller'
require_relative './server/controllers/upload_controller'
require_relative './server/controllers/download_controller'
require_relative './server/controllers/client_controller'
require_relative './server/controllers/user_log_controller'

class Metis
  class Server < Etna::Server
    get '/', action: 'client#index'
    get '/user', action: 'client#user'

    post '/check', action: 'user_log#check_log'

    post '/authorize/upload', action: 'upload#authorize', auth: { user: { can_edit?: :project_name } }

    post '/upload', action: 'upload#upload', auth: { hmac: true }

    #post '/upload/start', action: 'upload#start'
    #post '/upload/blob', action: 'upload#blob'
    #post '/upload/pause', action: 'upload#pause'
    #post '/upload/cancel', action: 'upload#cancel'
    #post '/upload/remove-file', action: 'upload#remove_file'
    #post '/upload/remove-failed', action: 'upload#remove_failed'
    #post '/upload/recover', action: 'upload#recover'

    get '/:project_name/download/*file_name', action: 'download#download', auth: { hmac: true }

    def initialize(config)
      super
      application.load_models
    end
  end
end
