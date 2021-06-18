require 'yaml'
require 'bundler'
Bundler.require(:default)

require_relative './lib/polyphemus'
require_relative './lib/server'

Polyphemus.instance.configure(YAML.load(File.read("config.yml")))

use Etna::CrossOrigin
use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes

run Polyphemus::Server.new
