require 'yaml'
require 'bundler'
Bundler.require(:default)

require './lib/metis'
require './lib/server'

use Rack::Static, urls: ['/css', '/fonts', '/img'], root: 'lib/client'
use Etna::ParseBody
use Etna::SymbolizeParams
use Metis::SetUid

use Etna::Auth
run Metis::Server.new(YAML.load(File.read("config.yml")))
