require_relative './server/controllers/vulcan_controller'
require_relative './server/controllers/browse_controller'
require_relative './server/controllers/workflows_controller'
require_relative './server/controllers/figure_controller'
require_relative './server/controllers/data_controller'
require_relative './server/controllers/sessions_controller'

class Vulcan
  class Server < Etna::Server
    # welcome_controller.rb
    get 'no_auth', as: :no_auth do
      erb_view(:no_auth)
    end

    get 'api/workflows', action: 'workflows#fetch', as: :workflows_view, auth: { user: { active?: true, has_flag?: 'vulcan' } }

    with auth: { user: { can_view?: :project_name, has_flag?: 'vulcan' } } do
      get 'api/:project_name/data/:cell_hash/:data_filename', action: 'data#fetch', as: :data_view, match_ext: true
      post 'api/:project_name/session/:workflow_name/status', action: 'sessions#status', as: :status_view, match_ext: true
      post 'api/:project_name/session/:workflow_name', action: 'sessions#submit', as: :submit_view, match_ext: true

      get 'api/:project_name/figures', action: 'figures#fetch'
      get 'api/:project_name/figure/:figure_id', action: 'figures#get'
      post 'api/:project_name/figure/:figure_id/update', action: 'figures#update'
    end

    with auth: { user: { active?: true, has_flag?: 'vulcan' } } do

      # remaining view routes are parsed by the client and must also be set there
      get '/:project_name/workflow', as: :workflow do
        erb_view(:client)
      end

      get '/:project_name/workflow/*view_path', as: :workflow_view do
        erb_view(:client)
      end
    end

    # root path
    get '/', as: :root do erb_view(:client) end
    get '/:project_name', as: :root do erb_view(:client) end

    def initialize
      super
      application.setup_db
    end
  end
end
