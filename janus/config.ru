# The packages
require 'yaml'
require 'bundler'
Bundler.require(:default, ENV["JANUS_ENV"].to_sym)

# The application
require_relative './lib/janus'
require_relative './lib/server'
require_relative './lib/server/throttle'
require_relative './lib/server/refresh_token'


Janus.instance.configure(YAML.load(File.read('config.yml')))

use Etna::CrossOrigin
use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'lib/client'
use Etna::Auth

use Janus::Throttle, max: 100
use Janus::RefreshToken
run Janus::Server.new
