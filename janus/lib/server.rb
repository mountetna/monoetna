# This class handles the http request and routing
require_relative './server/controllers/janus_controller'
require_relative './server/controllers/admin_controller'
require_relative './server/controllers/user_controller'
require_relative './server/controllers/authorization_controller'

class Janus
  class Server < Etna::Server
    get '/login', action: 'authorization#login', auth: { noauth: true }

    post '/api/validate-login', action: 'authorization#validate_login', auth: { noauth: true }

    # This generates nonces
    get '/time-signature', action: 'authorization#time_signature', auth: { noauth: true }

    # This uses a signed nonce to generate a new token
    get '/generate', action: 'authorization#generate', auth: { noauth: true }

    get '/api/tokens/nonce', action: 'authorization#time_signature', auth: { noauth: true }
    post '/api/tokens/generate', action: 'authorization#generate', auth: { noauth: true }
    post '/api/tokens/validate_task', action: 'authorization#validate_task', auth: { ignore_janus: true }

    get '/api/admin/:project_name/info', action: 'admin#project', auth: { user: { can_edit?: :project_name } }
    post '/api/admin/:project_name/update', action: 'admin#update_project', auth: { user: { is_supereditor?: true } }
    post '/api/admin/:project_name/update_permission', action: 'admin#update_permission', auth: { user: { is_admin?: :project_name } }
    post '/api/admin/:project_name/add_user', action: 'admin#add_user', auth: { user: { is_admin?: :project_name } }

    get '/api/admin/projects', action: 'admin#projects', auth: { user: { is_superviewer?: true } }
    post '/api/admin/add_project', action: 'admin#add_project', auth: { user: { is_supereditor?: true } }
    post '/api/admin/flag_user', action: 'admin#flag_user', auth: { user: { is_superuser?: true } }


    post '/api/user/update_key', action: 'user#update_key'
    get '/api/user/info', action: 'user#info'
    get '/api/user/projects', action: 'user#projects', auth: { user: { active?: true }, ignore_janus: true }
    get '/api/users', action: 'user#fetch_all', auth: { user: { is_superuser?: true } }

    get '/*views' do erb_view(:client) end

    get '/' do erb_view(:client) end

    def initialize
      super
      application.setup_db
    end
  end
end
