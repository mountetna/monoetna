require 'yaml'
require 'bundler'
Bundler.require(:default)

require './lib/metis'
require './lib/server'

Metis.instance.configure(YAML.load(File.read("config.yml")))

use Rack::Static, urls: ['/css', '/fonts', '/img'], root: 'lib/client'
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Metis::SetUid

run Metis::Server.new
