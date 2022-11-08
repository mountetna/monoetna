require 'yaml'
require 'bundler'
Bundler.require(:default)

require_relative './lib/gnomon'
require_relative './lib/server'

Gnomon.instance.configure(YAML.load(File.read("config.yml")))

use Etna::CrossOrigin
use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes

run Gnomon::Server.new
