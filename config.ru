# config.ru

require 'rack'
require 'json'
require 'digest'
require 'fileutils'
require 'redis'
require 'net/http'
require 'logger'

require './server/conf'
require './server/secrets'

require './server/controllers/controller'
require './server/controllers/client_controller'
require './server/controllers/log_controller'
require './server/controllers/upload_controller'
require './server/controllers/admin_controller'

require './server/service/sign_service'
require './server/service/redis_service'

require './server/metis'
require './server/routes'

# Set up access log
::Logger.class_eval { alias :write :'<<' } 
log_file = ::File.join(::File.dirname(::File.expand_path(__FILE__)),'log','access.log')
access_logger = ::Logger.new(log_file, 5, 1048576)
use ::Rack::CommonLogger, access_logger

use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'client'

run(Metis)