require 'json'
require 'rack'
require_relative '../magma'
require_relative '../magma/server/retrieve'
require_relative '../magma/server/query'
require_relative '../magma/server/update'
require_relative '../magma/server/gnomon'
require_relative '../magma/server/update_model'
require_relative '../magma/server/flag'

class Magma
  class Server < Etna::Server
    def initialize
      super
      application.tap do |magma|
        magma.load_models
      end
    end

    # Connect to the database and get some data.
    post '/retrieve', as: :retrieve, action: 'retrieve#action', auth: { user: { can_view?: :project_name } }

    post '/query', as: :query, action: 'query#action', auth: { user: { can_view?: :project_name } }

    post '/update', as: :update, action: 'update#action', auth: { user: { can_edit?: :project_name } }

    post '/update_model', action: 'update_model#action', auth: { user: { is_admin?: :project_name } }

    get '/flags/:project_name', action: 'flags#get', auth: { user: { can_view?: :project_name } }
    post '/flags/:project_name', action: 'flags#set', auth: { user: { is_admin?: :project_name } }

    post '/gnomon/rules', action: 'gnomon#rules', auth: { user: { is_supereditor?: true } }
    get '/gnomon/:project_name', action: 'gnomon#get', auth: { user: { can_view?: :project_name } }
    get '/gnomon/:project_name/revisions', action: 'gnomon#revisions', auth: { user: { can_view?: :project_name } }
    post '/gnomon/:project_name', action: 'gnomon#set', auth: { user: { is_admin?: :project_name } }
    post '/gnomon/:project_name/increment/:rule_name/:identifier_root', action: 'gnomon#increment', auth: { user: { is_admin?: :project_name } }
    get '/gnomon/:project_name/decompose/*identifier', action: 'gnomon#decompose', auth: { user: { can_view?: :project_name } }
    get '/gnomon/:project_name/list/:rule_name', action: 'gnomon#list', auth: { user: { can_view?: :project_name } }
    get '/gnomon/:project_name/rule/:rule_name', action: 'gnomon#rule', auth: { user: { can_view?: :project_name } }
    post '/gnomon/:project_name/generate/:rule_name/:identifier', action: 'gnomon#generate', auth: { user: { is_admin?: :project_name } }

    get '/' do
      [ 200, {}, [ 'Magma is available.' ] ]
    end
  end
end
