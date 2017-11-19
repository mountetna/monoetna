require 'bundler'
Bundler.require(:default)

require 'rack'
require 'json'
require 'pg'
require 'sequel'
require 'digest'
require 'fileutils'
require 'net/http'
require 'logger'
require 'openssl'
require 'securerandom' # Required by Phusion Passenger

# The details
require './lib/metis'
require './lib/server'
require './lib/server/conf'

# The application
require './lib/server/service/sign_service'
require './lib/server/controllers/metis_controller'
require './lib/server/controllers/main_controller'
require './lib/server/controllers/upload_controller'
require './lib/server/controllers/client_controller'
require './lib/server/controllers/user_log_controller'

use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'lib/client'
use Etna::ParseBody
use Etna::SymbolizeParams
run Metis::Server.new(YAML.load(File.read("config.yml")))
