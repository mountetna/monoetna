require 'bundler'
require 'rack'
Bundler.require(:default, :test)

ENV['POLYPHEMUS_ENV'] = 'test'

require 'webmock/rspec'

require 'simplecov'
SimpleCov.start

require 'yaml'
require_relative '../lib/server'
require_relative '../lib/polyphemus'

Polyphemus.instance.configure(YAML.load(File.read('config.yml')))

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Etna::DescribeRoutes
  run Polyphemus::Server.new
end

def auth_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
end

RSpec.configure do |config|
  config.mock_with :rspec do |mocks|
    mocks.verify_partial_doubles = true
  end

  config.shared_context_metadata_behavior = :apply_to_host_groups
  config.example_status_persistence_file_path = 'spec/examples.txt'
  #config.warnings = true
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post("/#{endpoint}", hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
end
