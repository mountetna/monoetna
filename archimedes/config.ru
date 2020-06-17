# This file is used by Rack-based servers to start the application.

require 'yaml'
require 'bundler'
Bundler.require(:default)

require_relative 'lib/archimedes'
require_relative 'lib/server'

use Etna::ParseBody
use Etna::SymbolizeParams
use Etna::Auth

run Archimedes::Server.new(YAML.load(File.read('config.yml')))
