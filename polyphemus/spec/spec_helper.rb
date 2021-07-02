require 'bundler'
require 'rack'
require 'rack/test'
Bundler.require(:default, :test)

ENV['POLYPHEMUS_ENV'] = 'test'

require 'webmock/rspec'
require 'database_cleaner'
require 'simplecov'
SimpleCov.start

require 'yaml'
require 'etna/spec/vcr'

require 'fileutils'
require 'timecop'

require_relative '../lib/server'
require_relative '../lib/polyphemus'

setup_base_vcr(__dir__)

Polyphemus.instance.configure(YAML.load(File.read('config.yml')))

METIS_HOST = Polyphemus.instance.config(:metis)[:host]
RELEASE_BUCKET = Polyphemus.instance.config(:metis)[:release_bucket]
RESTRICT_BUCKET = Polyphemus.instance.config(:metis)[:restrict_bucket]
MAGMA_HOST = Polyphemus.instance.config(:magma)[:host]
REDCAP_HOST = Polyphemus.instance.config(:redcap)[:host]
TEST_TOKEN = Polyphemus.instance.config(:polyphemus)[:token]

PROJECT = 'mvir1'
REDCAP_PROJECT_CONFIG_DIR = 'lib/etls/redcap/projects'

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Etna::DescribeRoutes
  run Polyphemus::Server.new
end

AUTH_USERS = {
  superuser: {
    email: 'zeus@twelve-labors.org', name: 'Zeus', perm: 'a:administration'
  },
  administrator: {
    email: 'hera@twelve-labors.org', name: 'Hera', perm: 'a:test', exp: 86401608136635
  },
  editor: {
    email: 'eurystheus@twelve-labors.org', name: 'Eurystheus', perm: 'e:labors'
  },
  privileged_editor: {
    email: 'copreus@twelve-labors.org', name: 'Copreus', perm: 'E:labors'
  },
  viewer: {
    email: 'hercules@twelve-labors.org', name: 'Hercules', perm: 'v:labors'
  },
  non_user: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: ''
  }
}

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

  config.order = :random

  config.before(:suite) do
    # DatabaseCleaner.strategy = :transaction
    DatabaseCleaner.clean_with(:truncation)
  end

  config.around(:each) do |example|
    # Unfortunately, DatabaseCleaner + Sequel does not properly handle the auto_savepointing, which means that
    # exceptions handled in rescue blocks do not behave correctly in tests (where as they would be fine outside of
    # tests).  Thus, we are forced to manually handle the transaction wrapping of examples manually to set this option.
    # See: http://sequel.jeremyevans.net/rdoc/files/doc/testing_rdoc.html#label-rspec+-3E-3D+2.8
    #      https://github.com/jeremyevans/sequel/issues/908#issuecomment-61217226
    Polyphemus.instance.db.transaction(:rollback=>:always, :auto_savepoint=>true){ example.run }
  end
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post("/#{endpoint}", hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
end

def stub_parent_exists(params={})
  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{params[:bucket] || RESTRICT_BUCKET}\//)
  .to_return({
    status: params[:status] || 200
  })
end

def stub_create_folder(params={})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{params[:bucket] || RESTRICT_BUCKET}\//)
  .to_return({
    status: params[:status] || 200
  })
end

def stub_rename_folder(params={})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{params[:bucket] || RESTRICT_BUCKET}\//)
  .to_return({
    status: params[:status] || 200
  })
end

def stub_magma_restricted_pools(base_model, restricted_pools)
  stub_request(:post, "#{MAGMA_HOST}/query")
    .with(body: hash_including({ project_name: 'mvir1',
                                query: [ base_model,
                                          [ 'timepoint', 'patient', 'restricted', '::true' ],
                                          '::all', "#{base_model}_pool", '::identifier' ] }))
    .to_return({ body: {
        'answer': restricted_pools.map {|p| [nil, p] }
    }.to_json })
end

def stub_magma_all_pools(base_model, all_pools)
  stub_request(:post, "#{MAGMA_HOST}/query")
    .with(body: hash_including({ project_name: 'mvir1',
                                query: [ "#{base_model}_pool", '::all', '::identifier' ] }))
    .to_return({ body: {
        'answer': all_pools.map {|p| [nil, p] }
    }.to_json })
end

def stub_magma_setup(patient_documents)
  stub_request(:post, "#{MAGMA_HOST}/retrieve")
    .with(body: hash_including({ project_name: 'mvir1', model_name: 'patient',
                                attribute_names: ['name', 'consent', 'restricted'], record_names: 'all' }))
    .to_return({ body: {
        'models': { 'patient': { 'documents': patient_documents } }
    }.to_json })

  stub_request(:post, "#{MAGMA_HOST}/update")
    .to_return do |request|

  body = StringIO.new(request.body)
  content_length = body.read.length
  body.rewind

  tempfile = Rack::Multipart::Parser::TEMPFILE_FACTORY
  bufsize = Rack::Multipart::Parser::BUFSIZE
  params = Rack::Utils.default_query_parser

  info = Rack::Multipart::Parser.parse body, content_length, request.headers['Content-Type'], tempfile, bufsize, params

  expect(info.params["project_name"]).to eq("mvir1")
  @all_updates << info.params["revisions"]
  { body: '{}' }
  end
end

def stub_metis_setup
  route_payload = JSON.generate([
    {:method=>"GET", :route=>"/:project_name/list_all_folders/:bucket_name", :name=>"folder_list_all_folders", :params=>["project_name", "bucket_name"]},
    {:method=>"GET", :route=>"/:project_name/list/:bucket_name/*folder_path", :name=>"folder_list", :params=>["project_name", "bucket_name", "folder_path"]},
    {:method=>"POST", :route=>"/:project_name/folder/rename/:bucket_name/*folder_path", :name=>"folder_rename", :params=>["project_name", "bucket_name", "folder_path"]},
    {:method=>"POST", :route=>"/:project_name/folder/create/:bucket_name/*folder_path", :name=>"folder_create", :params=>["project_name", "bucket_name", "folder_path"]}
  ])

  stub_request(:options, METIS_HOST).
    to_return({
      status: 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: route_payload
    })
  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    .to_return({
      status: 200,
      headers: {
      'Content-Type' => 'application/json'
      },
      body: JSON.parse(File.read('spec/fixtures/metis_release_folder_fixture.json')).to_json
    })

  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
    .to_return({
      status: 200,
      headers: {
      'Content-Type' => 'application/json'
      },
      body: JSON.parse(File.read('spec/fixtures/metis_restrict_folder_fixture.json')).to_json
    })
end

def stub_magma_models(fixture: 'spec/fixtures/magma_test_models.json')
  stub_request(:post, "#{MAGMA_HOST}/retrieve")
    .to_return({ body: File.read(fixture) })
end

def stub_magma_update_json
  stub_request(:post, "#{MAGMA_HOST}/update")
    .to_return({ status: 200, body: {}.to_json, headers: { 'Content-Type': 'application/json' } })
end

def stub_redcap_data(stub=nil)
  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: hash_including({ content: 'metadata' }))
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_metadata.json')
    })

  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: /fields/)
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_all.json')
    })
end

def stub_redcap_test2_data(stub=nil)
  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: hash_including({ content: 'metadata' }))
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_test2_metadata.json')
    })

  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: /fields/)
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_test2_data_all.json')
    })
end

def stub_redcap_multi_project_records
  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: /today/)
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_calendar.json')
    }).then
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_project_2_calendar.json')
    })

  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: /date_of_birth/)
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_essential_data.json')
    }).then
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_project_2_essential_data.json')
    })

  stub_request(:post, "#{REDCAP_HOST}/api/")
    .with(body: /height/)
    .to_return({
      headers: {
        'Content-Type': 'application/json'
      },
      body: File.read('spec/fixtures/redcap_mock_data_statistics.json')
    })
end

def copy_redcap_project(project_name='test')
  # Make sure the test project is in the right location so the
  # REDCap ETL can find it.
  test_fixture_path = "spec/fixtures/etls/redcap/#{project_name}.rb"
  test_output_path = "#{REDCAP_PROJECT_CONFIG_DIR}/#{project_name}.rb"
  FileUtils.mkdir_p(REDCAP_PROJECT_CONFIG_DIR) unless Dir.exist?(REDCAP_PROJECT_CONFIG_DIR)

  # Make sure we have the newest test project.
  File.delete(test_output_path) if File.file?(test_output_path)
  FileUtils.cp(test_fixture_path, test_output_path)
end

def temp_id(records, id)
  all_record_keys = []
  records.keys.each do |key|
    all_record_keys.concat(records[key].keys)
  end

  all_record_keys.each do |key|
    return key if key =~ /::temp-#{id}-.*/
  end
end
