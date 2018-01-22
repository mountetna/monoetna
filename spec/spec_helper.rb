require 'yaml'
require 'logger'
require 'rack/test'
require 'simplecov'
SimpleCov.start
require 'bundler'
Bundler.require(:default)

ENV['METIS_ENV'] = 'test'

require_relative '../lib/metis'
require_relative '../lib/server'

OUTER_APP = Rack::Builder.new do
  use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'lib/client'
  use Etna::ParseBody
  use Etna::SymbolizeParams

  use Etna::TestAuth
  run Metis::Server.new(YAML.load(File.read("config.yml")))
end

RSpec.configure do |config|
  config.expect_with :rspec do |expectations|
    # This option will default to `true` in RSpec 4. It makes the `description`
    # and `failure_message` of custom matchers include text for helper methods
    # defined using `chain`, e.g.:
    #     be_bigger_than(2).and_smaller_than(4).description
    #     # => "be bigger than 2 and smaller than 4"
    # ...rather than:
    #     # => "be bigger than 2"
    expectations.include_chain_clauses_in_custom_matcher_descriptions = true
  end

  config.mock_with :rspec do |mocks|
    # Prevents you from mocking or stubbing a method that does not exist on
    # a real object. This is generally recommended, and will default to
    # `true` in RSpec 4.
    mocks.verify_partial_doubles = true
  end

  # This option will default to `:apply_to_host_groups` in RSpec 4 (and will
  # have no way to turn it off -- the option exists only for backwards
  # compatibility in RSpec 3). It causes shared context metadata to be
  # inherited by the metadata hash of host groups and examples, rather than
  # triggering implicit auto-inclusion in groups with matching metadata.
  config.shared_context_metadata_behavior = :apply_to_host_groups

  # Allows RSpec to persist some state between runs in order to support
  # the `--only-failures` and `--next-failure` CLI options. We recommend
  # you configure your source control system to ignore this file.
  config.example_status_persistence_file_path = "spec/examples.txt"

  config.include FactoryBot::Syntax::Methods

  config.before(:suite) do
    FactoryBot.find_definitions

    DatabaseCleaner.strategy = :transaction
    DatabaseCleaner.clean_with(:truncation)
  end

  config.around(:each) do |example|
    DatabaseCleaner.cleaning do
      example.run
    end
  end
end

FactoryBot.define do
  factory :file, class: Metis::File do
    to_create(&:save)
  end
end

def fixture(name)
  File.join(File.dirname(__FILE__), "fixtures/#{name}.txt")
end

def json_body(body)
  JSON.parse(body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post("/#{endpoint}", hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
end
