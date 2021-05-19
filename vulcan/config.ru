# This file is used by Rack-based servers to start the application.

require 'yaml'
require 'bundler'
Bundler.require(:default)

require_relative 'lib/vulcan'
require_relative 'lib/server'

Vulcan.instance.configure(YAML.load(File.read('config.yml')))

use Etna::MetricsExporter
use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth
use Etna::DescribeRoutes

run Vulcan::Server.new
