require 'json'

require_relative './server/controllers/metis_controller'
require_relative './server/controllers/files_controller'
require_relative './server/controllers/upload_controller'
require_relative './server/controllers/download_controller'
require_relative './server/controllers/client_controller'
require_relative './server/set_uid'

class Metis
  class Server < Etna::Server
    get '/:project_name', action: 'client#index', auth: { user: { can_view?: :project_name } }

    post '/authorize/upload', action: 'upload#authorize', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/upload/*file_name', action: 'upload#upload', auth: { hmac: true }, as: :upload

    get '/:project_name/download/*file_name', action: 'download#download', auth: { hmac: true }

    get '/:project_name/files', action: 'files#index', auth: { user: { can_view?: :project_name } }

    def initialize(config)
      super
      application.load_models
    end
  end
end
