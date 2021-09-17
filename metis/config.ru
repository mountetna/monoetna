require 'yaml'
require 'bundler'
Bundler.require(:default)

require './lib/metis'
require './lib/server'

Metis.instance.configure(YAML.load(File.read("config.yml")))

use Etna::CrossOrigin
use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes
use Metis::SetUid
use Etna::SynchronizeDb

run Metis::Server.new
