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
    email: 'hera@olympus.org', name: 'Hera', perm: 'a:labors', exp: Time.now.to_i + 6000, flags: 'vulcan'
  },
  editor: {
    email: 'eurystheus@twelve-labors.org', name: 'Eurystheus', perm: 'e:labors', exp: Time.now.to_i + 6000, flags: 'vulcan'
  },
  viewer: {
    email: 'hercules@twelve-labors.org', name: 'Hercules', perm: 'v:labors', exp: Time.now.to_i + 6000, flags: 'vulcan'
  },
  non_user: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: '', exp: Time.now.to_i + 6000, flags: 'vulcan'
  },
  no_flag: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: 'g:labors', exp: Time.now.to_i + 6000,
  },
  guest: {
    email: 'sinon@troy.org', name: 'Sinon', perm: 'g:labors', exp: Time.now.to_i + 6000, flags: 'vulcan'
  }
}

PROJECT = "labors"

def store(hash, filename, data)
  storage = Vulcan::Storage.new

  path = storage.data_path(project_name: PROJECT, cell_hash: hash, data_filename: filename)
  ::FileUtils.mkdir_p(::File.dirname(path))
  ::File.write(path, data)
  path
end

def clear_store
  storage = Vulcan::Storage.new

  FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)
end


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

def create_figure(params)
  now = DateTime.now
  create(
    :figure,
    {
      figure_id: 1,
      project_name: 'labors',
      workflow_name: 'workflow',
      author: 'author',
      inputs: {},
      title: 'title',
      created_at: now,
      updated_at: now
    }.update(params)
  )
end

FactoryBot.define do
  factory :figure, class: Vulcan::Figure do
    to_create(&:save)
  end

  factory :workflow_snapshot, class: Vulcan::WorkflowSnapshot do
    to_create(&:save)
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

def create_figure_with_snapshot
  fig = create_figure(
    title: "Lion of Nemea",
    workflow_name: "test_workflow.cwl",
    dependencies: {
      something: "sha:abc",
    },
  )
  
  # Tweak the snapshot to be slightly different
  fig.workflow_snapshot.update(
    inputs: {
      "someInt" => {
        "default" => 200,
        "format" => nil,
        "label" => "it is an int",
        "type" => "int",
        "doc" => "help tip",
      },
      "someIntWithoutDefault" => {
        "default" => nil,
        "format" => nil,
        "label" => nil,
        "type" => "int",
        "doc" => "another tip",
      },
      "removedInt" => {
        "default" => nil,
        "format" => nil,
        "label" => nil,
        "type" => "int",
        "doc" => "deprecated tip",
      },
    }, steps: [
      [
        {
          "in" => [{ "id" => "a", "source" => "someInt" },
                  { "id" => "b", "source" => "someIntWithoutDefault" }],
          "doc" => nil,
          "label" => nil,
          "out" => ["sum"],
          "id" => "firstAdd",
          "run" => "scripts/add.cwl",
        },
        {
          "in" => [{ "id" => "a", "source" => "firstAdd/sum" },
                  { "id" => "b", "source" => "removedInt" }],
          "doc" => nil,
          "label" => nil,
          "out" => ["sum"],
          "id" => "secondAdd",
          "run" => "scripts/add.cwl",
        },
        {
          "in" => [{ "id" => "num", "source" => "secondAdd/sum" }],
          "doc" => nil,
          "label" => nil,
          "out" => ["num"],
          "id" => "pickANum",
          "run" => "ui-queries/pick-a-number.cwl",
        },
        {
          "in" => [{ "id" => "a", "source" => "secondAdd/sum" },
                  { "id" => "b", "source" => "pickANum/num" }],
          "doc" => nil,
          "label" => nil,
          "out" => ["sum", "thumb.png"],
          "id" => "finalStep",
          "run" => "scripts/add.cwl",
        },
        {
          "in" => [{ "id" => "a", "source" => "finalStep/sum" },
                  { "id" => "b", "source" => "finalStep/thumb.png" }],
          "doc" => nil,
          "label" => nil,
          "id" => "aPlot",
          "out" => [],
          "run" => "ui-outputs/plotly.cwl",
        },
      ],
    ],
  )

  fig.refresh

  fig
end

def below_admin_roles
  [:editor, :viewer, :guest]
end

def below_editor_roles
  [:viewer, :guest]
end