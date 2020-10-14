require 'yaml'
require 'bundler'
Bundler.require(:default)

require './lib/polyphemus'
require './lib/server'

Polyphemus.instance.configure(YAML.load(File.read("config.yml")))

use Etna::CrossOrigin
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes

run Polyphemus::Server.new
