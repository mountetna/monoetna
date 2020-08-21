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

METIS_HOST = Polyphemus.instance.config(:metis)[:host]
RELEASE_BUCKET = Polyphemus.instance.config(:metis)[:release_bucket]
RESTRICT_BUCKET = Polyphemus.instance.config(:metis)[:restrict_bucket]
PROJECT = 'mvir1'

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

def stub_magma_linked_pools(patient_name, linked_pools)
  stub_request(:post, 'https://magma.test/query')
    .with(body: hash_including({ project_name: 'mvir1',
                                query: [ 'cytof',
                                          [ 'timepoint', 'patient', 'name', '::equals', patient_name ],
                                          '::all', 'cytof_pool', '::identifier' ] }))
    .to_return({ body: {
        'answer': linked_pools.map {|p| [nil, p] }
    }.to_json })
end

def stub_magma_pool_state(pool_name, pool_states)
  stub_request(:post, 'https://magma.test/query')
    .with(body: hash_including({ project_name: 'mvir1',
                                query: [ 'cytof_pool', ['pool_name', '::equals', pool_name], '::all', 'restricted' ] }))
    .to_return({ body: {
        'answer': pool_states.map {|p| [nil, p] }
    }.to_json })
end

def stub_magma_setup(patient_documents)
  stub_request(:post, 'https://magma.test/retrieve')
    .with(body: hash_including({ project_name: 'mvir1', model_name: 'patient',
                                attribute_names: 'all', record_names: 'all' }))
    .to_return({ body: {
        'models': { 'patient': { 'documents': patient_documents } }
    }.to_json })

  stub_request(:post, 'https://magma.test/update')
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