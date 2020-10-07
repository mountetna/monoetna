require 'yaml'
require 'logger'
require 'simplecov'
require 'rack/test'
require 'bundler'
require 'securerandom'
require 'timecop'
require 'webmock/rspec'
require 'vcr'

Bundler.setup(:default, :test)

SimpleCov.start

ENV['ARACHNE_ENV'] = 'test'

require_relative '../lib/etna'

def setup_app(server, layer=nil, config={ test: {} })
  Etna::Application.find(server).configure(config)
  Rack::Builder.new do
    use Etna::ParseBody
    use Etna::SymbolizeParams
    use *layer if layer
    run server.new
  end
end

VCR.configure do |c|
  c.hook_into :webmock
  c.cassette_serializers
  c.cassette_library_dir = ::File.join(__dir__,  'fixtures', 'cassettes')
  c.allow_http_connections_when_no_cassette = true
  c.register_request_matcher :try_body do |request_1, request_2|
    if request_1.headers['Content-Type'].first =~ /application\/json/
      if request_2.headers['Content-Type'].first =~ /application\/json/
        JSON.parse(request_1.body) == JSON.parse(request_2.body)
      else
        false
      end
    else
      request_1.body == request_2.body
    end
  end

  c.default_cassette_options = {
      serialize_with: :compressed,
      record: if ENV['IS_CI']
        :none
      else
        ENV['RERECORD'] ? :all : :once
      end,
      match_requests_on: [:method, :uri, :try_body]
  }

  # Filter the authorization headers of any request by replacing any occurrence of that request's
  # Authorization value with <AUTHORIZATION>
  c.filter_sensitive_data('<AUTHORIZATION>') do |interaction|
    interaction.request.headers['Authorization'].first
  end

  # The model synchronization workflow benefits from testing against production models, but produces a massive
  # recording file that doesn't fit into git easily.  To assist in that, we prune some large chunks of data
  # that are generally not useful to recording tests.  But this is fragile and assumes that no tests using vcr depend
  # on certain response contents.  Ideal solution would be to greatly reduce the size of default magma retrieve and
  # update_model responses to not be as noisy for large models.
  c.filter_sensitive_data('{}') do |interaction|
    if interaction.request.uri =~ /\/update_model/
      interaction.response.body
    else
      "{}"
    end
  end

  c.filter_sensitive_data('"options":["one-option"]') do |interaction|
    match = interaction.response.body.match(/"options":\["ENSG0[^\]]+\]/)
    match ? match[0] : '"options":["one-option"]'
  end

  c.filter_sensitive_data('"value":["one-option"]') do |interaction|
    match = interaction.response.body.match(/"value":\["ENSG0[^\]]+\]/)
    match ? match[0] : '"value":["one-option"]'
  end
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

  # rspec-mocks config goes here. You can use an alternate test double
  # library (such as bogus or mocha) by changing the `mock_with` option here.
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
end

METIS_HOST = 'https://metis.test'
JANUS_HOST = 'https://janus.test'
MAGMA_HOST = 'https://magma.test'
PROJECT = 'test'
RESTRICT_BUCKET = 'restrict'
RELEASE_BUCKET = 'release'

def stub_metis_setup
  route_payload = JSON.generate([
    {:method=>"GET", :route=>"/:project_name/list_all_folders/:bucket_name", :name=>"folder_list_all_folders", :params=>["project_name", "bucket_name"]},
    {:method=>"GET", :route=>"/:project_name/list/:bucket_name/*folder_path", :name=>"folder_list", :params=>["project_name", "bucket_name", "folder_path"]},
    {:method=>"POST", :route=>"/:project_name/folder/rename/:bucket_name/*folder_path", :name=>"folder_rename", :params=>["project_name", "bucket_name", "folder_path"]},
    {:method=>"POST", :route=>"/:project_name/folder/create/:bucket_name/*folder_path", :name=>"folder_create", :params=>["project_name", "bucket_name", "folder_path"]},
    {:method=>"POST", :route=>"/:project_name/find/:bucket_name", :name=>"bucket_find", :params=>["project_name", "bucket_name"]},
    {:method=>"POST", :route=>"/:project_name/files/copy", :name=>"file_bulk_copy", :params=>["project_name"]}
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

def stub_find(params={})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/find\/#{params[:bucket] || RESTRICT_BUCKET}/)
  .to_return({
    status: params[:status] || 200
  })
end

def stub_copy(params={})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/files\/copy/)
  .to_return({
    status: params[:status] || 200
  })
end

def stub_janus_setup
  stub_request(:get, /#{JANUS_HOST}\/refresh_token/)
    .to_return({
      status: 200,
      body: 'a token for you!'
    })

  stub_request(:get, /#{JANUS_HOST}\/project\/#{PROJECT}/)
    .to_return({
      status: 200,
      body: '<html><body>A project</body></html>'
    })

  stub_request(:get, /#{JANUS_HOST}\/whoami/)
    .to_return({
      status: 200,
      body: '{"email": "janus@twofaces.org"}'
    })

  stub_request(:post, /#{JANUS_HOST}\/add_project/)
    .to_return({
      status: 302
    })

  stub_request(:post, /#{JANUS_HOST}\/add_user/)
    .to_return({
      status: 302
    })

  stub_request(:post, /#{JANUS_HOST}\/update_permission/)
    .to_return({
      status: 302
    })

  stub_request(:get, /#{JANUS_HOST}\/viewer_token/)
    .to_return({
      status: 200,
      body: 'a view-only token for you!'
    })
end

def stub_magma_models(models)
  stub_request(:post, /#{MAGMA_HOST}\/retrieve/)
  .to_return({
    headers: {
      'Content-Type': 'application/json'
    },
    body: models.to_json
  })
end

def stub_magma_update
  stub_request(:post, /#{MAGMA_HOST}\/update/)
  .to_return do |request|

    body = StringIO.new(request.body)
    content_length = body.read.length
    body.rewind

    tempfile = Rack::Multipart::Parser::TEMPFILE_FACTORY
    bufsize = Rack::Multipart::Parser::BUFSIZE
    params = Rack::Utils.default_query_parser

    info = Rack::Multipart::Parser.parse body, content_length, request.headers['Content-Type'], tempfile, bufsize, params

    expect(info.params["project_name"]).to eq(PROJECT)
    @all_updates << info.params["revisions"]
    { body: '{}' }
    end
end

def stub_magma_update_model
  stub_request(:post, /#{MAGMA_HOST}\/update_model/)
  .to_return({
    status: 200,
    body: {}.to_json
  })
end