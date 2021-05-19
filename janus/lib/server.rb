# This class handles the http request and routing
require_relative './server/controllers/janus_controller'
require_relative './server/controllers/admin_controller'
require_relative './server/controllers/user_controller'
require_relative './server/controllers/authorization_controller'

class Janus
  class Server < Etna::Server
    get '/login', action: 'authorization#login', auth: { noauth: true }

    post '/validate-login', action: 'authorization#validate_login', auth: { noauth: true }

    # This generates nonces
    get '/time-signature', action: 'authorization#time_signature', auth: { noauth: true }

    # This uses a signed nonce to generate a new token
    get '/generate', action: 'authorization#generate', auth: { noauth: true }

    get '/refresh_token', action: 'user#refresh_token', auth: { user: { active?: true } }

    get '/api/tokens/nonce', action: 'authorization#time_signature', auth: { noauth: true }
    post '/api/tokens/generate', action: 'authorization#generate', auth: { noauth: true }
    post '/api/tokens/validate_task', action: 'authorization#validate_task', auth: { ignore_janus: true }

    # Once we figure out a long-term token strategy, this should probably get
    #   consolidated with /refresh_token so we don't just keep creating
    #   small token-related views.
    get '/viewer_token', action: 'user#viewer_token', auth: { user: { is_superuser?: true } }

    post '/update_permission/:project_name', action: 'admin#update_permission', auth: { user: { is_admin?: :project_name } }

    post '/add_user/:project_name', action: 'admin#add_user', auth: { user: { is_admin?: :project_name } }

    post '/add_project', action: 'admin#add_project', auth: { user: { is_superuser?: true } }

    post '/flag_user', action: 'admin#flag_user', auth: { user: { is_superuser?: true } }

    post '/update_key', action: 'user#update_key'

    get '/user', action: 'user#info'

    get '/users', action: 'user#fetch_all', auth: { user: { is_superuser?: true } }

    get '/allprojects', action: 'admin#projects', auth: { user: { is_superviewer?: true } }

    get '/projects', action: 'user#projects', auth: { user: { active?: true } }

    get '/project/:project_name', action: 'admin#project', auth: { user: { can_edit?: :project_name } }

    get '/admin' do erb_view(:client) end

    get '/settings' do erb_view(:client) end

    get '/:project_name', auth: { user: { can_edit?: :project_name } } do erb_view(:client) end

    get '/' do erb_view(:client) end

    def initialize
      super
      application.setup_db
    end
  end
end
