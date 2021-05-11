require 'bundler'
Bundler.require(:default, :test)

ENV['VULCAN_ENV'] = 'test'

require 'webmock/rspec'

require 'simplecov'
SimpleCov.start

require 'yaml'
require 'factory_bot'
require 'database_cleaner'
require 'rack/test'

require_relative '../lib/server'
require_relative '../lib/vulcan'
require 'etna/spec/vcr'

Vulcan.instance.configure(YAML.load(File.read('config.yml')))
setup_base_vcr(__dir__)

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Etna::DescribeRoutes
  run Vulcan::Server.new
end

AUTH_USERS = {
  admin: {
    email: 'hera@olympus.org', name: 'Hera', perm: 'a:labors', exp: Time.now.to_i + 6000,
  },
  editor: {
    email: 'eurystheus@twelve-labors.org', name: 'Eurystheus', perm: 'e:labors', exp: Time.now.to_i + 6000,
  },
  viewer: {
    email: 'hercules@twelve-labors.org', name: 'Hercules', perm: 'v:labors', exp: Time.now.to_i + 6000,
  },
  non_user: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: '', exp: Time.now.to_i + 6000,
  }
}

PROJECT = "labors"

def auth_header(user_type, task: false, additional: {})
  user = AUTH_USERS[user_type].dup
  user[:task] = task if task
  header(*Etna::TestAuth.token_header({}.update(user).update(additional)))
end

RSpec.configure do |config|
  config.mock_with :rspec do |mocks|
    mocks.verify_partial_doubles = true
  end

  config.shared_context_metadata_behavior = :apply_to_host_groups
  config.example_status_persistence_file_path = 'spec/examples.txt'
  #config.warnings = true

  config.include FactoryBot::Syntax::Methods

  config.before(:suite) do
    FactoryBot.find_definitions
    # DatabaseCleaner.strategy = :transaction
    DatabaseCleaner.clean_with(:truncation)
  end

  config.around(:each) do |example|
    # Unfortunately, DatabaseCleaner + Sequel does not properly handle the auto_savepointing, which means that
    # exceptions handled in rescue blocks do not behave correctly in tests (where as they would be fine outside of
    # tests).  Thus, we are forced to manually handle the transaction wrapping of examples manually to set this option.
    # See: http://sequel.jeremyevans.net/rdoc/files/doc/testing_rdoc.html#label-rspec+-3E-3D+2.8
    #      https://github.com/jeremyevans/sequel/issues/908#issuecomment-61217226
    Vulcan.instance.db.transaction(:rollback=>:always, :auto_savepoint=>true){ example.run }
  end
end


def run_script script
  txt = script

  manifest = Archimedes::Manifest.new(
    'xyzzy',
    'timur',
    txt
  )
  manifest.payload
  manifest.instance_variable_get('@return_vars')
end

FactoryBot.define do
  factory :view do
    to_create(&:save)
  end

  factory :manifest do
    to_create(&:save)
    project { 'labors' }
    sequence :name do |n|
      "manifest #{n}"
    end

    trait :script do
      script { '@value = 1 + 1' }
    end

    trait :public do
      access { 'public' }
    end

    trait :private do
      access { 'private' }
    end
  end

  factory :plot do
    to_create(&:save)
    project { 'labors' }
    script { '@test = 1' }
    sequence :name do |n|
      "plot #{n}"
    end

    trait :scatter do
      plot_type { 'scatter' }
      configuration { {plot: 'ok'} }
    end

    trait :public do
      access { 'public' }
    end

    trait :private do
      access { 'private' }
    end
  end

  factory :user do
    to_create(&:save)

    AUTH_USERS.each do |user_type, template|
      trait user_type do
        email { template[:email] }
        name { template[:name] }
      end
    end
  end
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post("/#{endpoint}", hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
end

def save_last_response_json(fixture_name, type)
  fixture_path = ::File.join(__dir__, '..', 'lib', 'client', 'jsx', 'test_utils', 'fixtures', "#{fixture_name}.ts")
  constName = fixture_name.gsub(/[_-](\w)/){$1.upcase}
  p last_response.body

  ::File.write(fixture_path, "import {#{type}} from \"../../api_types\";\n\nexport const #{constName}: #{type} = #{last_response.body};")
end