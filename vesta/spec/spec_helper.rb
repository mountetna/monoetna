require "bundler"
require "rack"
require "rack/test"
Bundler.require(:default, :test)

ENV["VESTA_ENV"] = "test"

require "webmock/rspec"
require "database_cleaner"
require "factory_bot"
require "simplecov"
SimpleCov.start

require "yaml"
require "etna/spec/vcr"

require "fileutils"
require "timecop"

require_relative "../lib/server"
require_relative "../lib/vesta"

Vesta.instance.configure(YAML.load(File.read("config.yml")))

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Etna::DescribeRoutes
  run Vesta::Server.new
end

RSpec.configure do |config|
  config.mock_with :rspec do |mocks|
    mocks.verify_partial_doubles = true
  end

  config.expect_with :rspec do |exp|
    exp.max_formatted_output_length = nil
  end

  config.shared_context_metadata_behavior = :apply_to_host_groups
  config.example_status_persistence_file_path = "spec/examples.txt"
  #config.warnings = true

  config.order = :random

  config.include FactoryBot::Syntax::Methods

  config.before(:suite) do
    FactoryBot.find_definitions
  end

  config.around(:each) do |example|
    # Unfortunately, DatabaseCleaner + Sequel does not properly handle the auto_savepointing, which means that
    # exceptions handled in rescue blocks do not behave correctly in tests (where as they would be fine outside of
    # tests).  Thus, we are forced to manually handle the transaction wrapping of examples manually to set this option.
    # See: http://sequel.jeremyevans.net/rdoc/files/doc/testing_rdoc.html#label-rspec+-3E-3D+2.8
    #      https://github.com/jeremyevans/sequel/issues/908#issuecomment-61217226
    Vesta.instance.db.transaction(:rollback => :always, :auto_savepoint => true) { example.run }
  end

  config.after(:each) do 
  end
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post(URI.encode("/#{endpoint.reverse.chomp('/').reverse}"), hash.to_json, { "CONTENT_TYPE" => "application/json" })
end

