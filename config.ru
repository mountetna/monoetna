# config.ru

require 'rack'
require 'json'
require 'digest'
require 'fileutils'
require 'redis'
require 'net/http'
require 'logger'

# The details
require './server/conf'
require './server/secrets'

# The database
require './server/service/sign_service'
require './server/service/redis_service'

# The application
require './server/errors/basic_error'
require './server/models/user_model'
require './server/models/file_model'
require './server/metis'
require './server/routes'
require './server/controllers/basic_controller'
require './server/controllers/main_controller'
require './server/controllers/upload_controller'
require './server/controllers/client_controller'
require './server/controllers/user_log_controller'

#require './server/controllers/controller'
#require './server/controllers/upload_controller'

use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'client'
run(Metis)