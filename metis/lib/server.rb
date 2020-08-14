require 'json'

require_relative './server/controllers/metis_controller'
require_relative './server/controllers/file_controller'
require_relative './server/controllers/folder_controller'
require_relative './server/controllers/bucket_controller'
require_relative './server/controllers/upload_controller'
require_relative './server/controllers/download_controller'
require_relative './server/controllers/client_controller'
require_relative './server/set_uid'

class Metis
  class Server < Etna::Server
    get '/', action: 'client#index'
    get '/:project_name', action: 'client#index', auth: { user: { can_view?: :project_name } }
    get '/:project_name/browse/:bucket_name', action: 'client#index', auth: { user: { can_view?: :project_name } }
    get '/:project_name/browse/:bucket_name/*folder_name', action: 'client#index', auth: { user: { can_view?: :project_name } }

    post '/authorize/upload', action: 'upload#authorize', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/upload/:bucket_name/*file_path', action: 'upload#upload', auth: { hmac: true }, as: :upload

    post '/authorize/download', action: 'download#authorize', auth: { user: { can_view?: :project_name } }
    get '/:project_name/download/:bucket_name/*file_path', action: 'download#download', auth: { hmac: true }, as: :download

    get '/:project_name/list/:bucket_name/*folder_path', action: 'folder#list', auth: { user: { can_view?: :project_name } }
    get '/:project_name/list/:bucket_name', action: 'folder#list', auth: { user: { can_view?: :project_name } }
    get '/:project_name/list/', action: 'bucket#list', auth: { user: { can_view?: :project_name } }

    # bucket operations
    post '/:project_name/bucket/create/:bucket_name', action: 'bucket#create', auth: { user: { is_admin?: :project_name } }
    post '/:project_name/bucket/update/:bucket_name', action: 'bucket#update', auth: { user: { is_admin?: :project_name } }
    delete '/:project_name/bucket/remove/:bucket_name', action: 'bucket#remove', auth: { user: { is_admin?: :project_name } }

    # folder operations
    post '/:project_name/folder/create/:bucket_name/*folder_path', action: 'folder#create', auth: { user: { can_edit?: :project_name } }
    delete '/:project_name/folder/remove/:bucket_name/*folder_path', action: 'folder#remove', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/folder/protect/:bucket_name/*folder_path', action: 'folder#protect', auth: { user: { is_admin?: :project_name } }
    post '/:project_name/folder/unprotect/:bucket_name/*folder_path', action: 'folder#unprotect', auth: { user: { is_admin?: :project_name } }
    post '/:project_name/folder/rename/:bucket_name/*folder_path', action: 'folder#rename', auth: { user: { can_edit?: :project_name } }

    # file operations
    delete '/:project_name/file/remove/:bucket_name/*file_path', action: 'file#remove', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/file/protect/:bucket_name/*file_path', action: 'file#protect', auth: { user: { is_admin?: :project_name } }
    post '/:project_name/file/unprotect/:bucket_name/*file_path', action: 'file#unprotect', auth: { user: { is_admin?: :project_name } }
    post '/:project_name/file/rename/:bucket_name/*file_path', action: 'file#rename', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/file/copy/:bucket_name/*file_path', action: 'file#copy', auth: { user: { can_edit?: :project_name } }
    post '/:project_name/files/copy', action: 'file#bulk_copy', auth: { user: { can_edit?: :project_name } }

    def initialize
      super
      application.load_models
    end
  end
end
