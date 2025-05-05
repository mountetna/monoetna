require 'yaml'
require 'bundler'
Bundler.require(:default)

require_relative './lib/vesta'
require_relative './lib/server'

Vesta.instance.configure(YAML.load(File.read("config.yml")))

use Etna::CrossOrigin
use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes
use Etna::SynchronizeDb

run Vesta::Server.new
