# config.ru

require 'rack'
require 'json'
require 'digest'
require 'fileutils'
require 'redis'

require './server/conf'

require './server/controllers/controller'
require './server/controllers/client_controller'
require './server/controllers/upload_controller'

require './server/service/sign_service'
require './server/service/redis_service'

require './server/metis'
require './server/routes'

use Rack::Static, urls: ['/css', '/js', '/fonts'], root: 'client'

run(Metis)