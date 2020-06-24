require 'yaml'
require 'bundler'
Bundler.require(:default)

require './lib/metis'
require './lib/server'

Metis.instance.configure(YAML.load(File.read("config.yml")))

# Used primarily by development servers.  These are served via apache in production.
use Rack::Static, urls: ['/css', '/fonts', '/img', '/js'], root: 'public'
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes
use Metis::SetUid

run Metis::Server.new
