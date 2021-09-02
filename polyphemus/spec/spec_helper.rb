require "bundler"
require "rack"
require "rack/test"
Bundler.require(:default, :test)

ENV["POLYPHEMUS_ENV"] = "test"

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
require_relative "../lib/polyphemus"

setup_base_vcr(__dir__)

Polyphemus.instance.configure(YAML.load(File.read("config.yml")))

METIS_HOST = Polyphemus.instance.config(:metis)[:host]
RELEASE_BUCKET = Polyphemus.instance.config(:metis)[:release_bucket]
RESTRICT_BUCKET = Polyphemus.instance.config(:metis)[:restrict_bucket]
MAGMA_HOST = Polyphemus.instance.config(:magma)[:host]
REDCAP_HOST = Polyphemus.instance.config(:redcap)[:host]
TEST_TOKEN = Polyphemus.instance.config(:polyphemus)[:token]

PROJECT = "mvir1"
REDCAP_PROJECT_CONFIG_DIR = "lib/etls/redcap/projects"
RENAMING_PROJECT_CONFIG_DIR = "lib/etls/renaming/projects"

OUTER_APP = Rack::Builder.new do
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Etna::DescribeRoutes
  run Polyphemus::Server.new
end

AUTH_USERS = {
  superuser: {
    email: "zeus@twelve-labors.org", name: "Zeus", perm: "a:administration",
  },
  administrator: {
    email: "hera@twelve-labors.org", name: "Hera", perm: "a:test", exp: 86401608136635,
  },
  editor: {
    email: "eurystheus@twelve-labors.org", name: "Eurystheus", perm: "e:labors",
  },
  privileged_editor: {
    email: "copreus@twelve-labors.org", name: "Copreus", perm: "E:labors",
  },
  viewer: {
    email: "hercules@twelve-labors.org", name: "Hercules", perm: "v:labors",
  },
  non_user: {
    email: "nessus@centaurs.org", name: "Nessus", perm: "",
  },
}

def auth_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
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
    # DatabaseCleaner.strategy = :transaction
    DatabaseCleaner.clean_with(:truncation)
  end

  config.around(:each) do |example|
    # Unfortunately, DatabaseCleaner + Sequel does not properly handle the auto_savepointing, which means that
    # exceptions handled in rescue blocks do not behave correctly in tests (where as they would be fine outside of
    # tests).  Thus, we are forced to manually handle the transaction wrapping of examples manually to set this option.
    # See: http://sequel.jeremyevans.net/rdoc/files/doc/testing_rdoc.html#label-rspec+-3E-3D+2.8
    #      https://github.com/jeremyevans/sequel/issues/908#issuecomment-61217226
    Polyphemus.instance.db.transaction(:rollback => :always, :auto_savepoint => true) { example.run }
  end
end

FactoryBot.define do
  factory :ingest_file, class: Polyphemus::IngestFile do
    to_create(&:save)
  end

  factory :watch_folder, class: Polyphemus::WatchFolder do
    to_create(&:save)
  end
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post("/#{endpoint}", hash.to_json, { "CONTENT_TYPE" => "application/json" })
end

def stub_parent_exists(params = {})
  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
    })
end

def stub_create_folder(params = {})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
    })
end

def stub_rename_folder(params = {})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
    })
end

def stub_rename_folder_with_error(params = {})
  stub_request(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: 422,
    }).then.to_return({
      status: params[:status] || 200,
    })
end

def stub_magma_restricted_pools(base_model, restricted_pools)
  stub_request(:post, "#{MAGMA_HOST}/query")
    .with(body: hash_including({ project_name: "mvir1",
                                query: [base_model,
                                        ["timepoint", "patient", "restricted", "::true"],
                                        "::all", "#{base_model}_pool", "::identifier"] }))
    .to_return({ body: {
      'answer': restricted_pools.map { |p| [nil, p] },
    }.to_json })
end

def stub_magma_all_pools(base_model, all_pools)
  stub_request(:post, "#{MAGMA_HOST}/query")
    .with(body: hash_including({ project_name: "mvir1",
                                 query: ["#{base_model}_pool", "::all", "::identifier"] }))
    .to_return({ body: {
      'answer': all_pools.map { |p| [nil, p] },
    }.to_json })
end

def stub_magma_setup(patient_documents)
  stub_request(:post, "#{MAGMA_HOST}/retrieve")
    .with(body: hash_including({ project_name: "mvir1", model_name: "patient",
                                 attribute_names: ["name", "consent", "restricted"], record_names: "all" }))
    .to_return({ body: {
      'models': { 'patient': { 'documents': patient_documents } },
    }.to_json })

  stub_magma_update("mvir1")
end

def stub_magma_update(project_name = nil)
  stub_request(:post, /#{MAGMA_HOST}\/update$/)
    .to_return do |request|
    body = StringIO.new(request.body)
    content_length = body.read.length
    body.rewind

    tempfile = Rack::Multipart::Parser::TEMPFILE_FACTORY
    bufsize = Rack::Multipart::Parser::BUFSIZE
    params = Rack::Utils.default_query_parser

    info = Rack::Multipart::Parser.parse body, content_length, request.headers["Content-Type"], tempfile, bufsize, params

    expect(info.params["project_name"]).to eq(project_name) if project_name
    @all_updates << info.params["revisions"] if info.params && info.params["revisions"]
    { body: "{}" }
  end
end

def stub_metis_setup
  route_payload = JSON.generate([
    { :method => "GET", :route => "/:project_name/list_all_folders/:bucket_name", :name => "folder_list_all_folders", :params => ["project_name", "bucket_name"] },
    { :method => "GET", :route => "/:project_name/list/:bucket_name/*folder_path", :name => "folder_list", :params => ["project_name", "bucket_name", "folder_path"] },
    { :method => "GET", :route => "/:project_name/list_by_id/:bucket_name/:folder_id", :name => "folder_list_by_id", :params => ["project_name", "bucket_name", "folder_id"] },
    { :method => "GET", :route => "/:project_name/touch/:bucket_name/*folder_path", :name => "folder_touch", :params => ["project_name", "bucket_name", "folder_path"] },
    { :method => "POST", :route => "/:project_name/folder/rename/:bucket_name/*folder_path", :name => "folder_rename", :params => ["project_name", "bucket_name", "folder_path"] },
    { :method => "POST", :route => "/:project_name/folder/create/:bucket_name/*folder_path", :name => "folder_create", :params => ["project_name", "bucket_name", "folder_path"] },
    { :method => "POST", :route => "/authorize/upload", :name => "upload_authorize", :params => ["project_name", "bucket_name", "file_path"] },
    { :method => "POST", :route => "/:project_name/upload/:bucket_name/*file_path", :name => "upload_upload", :params => ["project_name", "bucket_name", "file_path"] },
    { :method => "POST", :route => "/:project_name/find/:bucket_name", :name => "bucket_find", :params => ["project_name", "bucket_name"] },
  ])

  stub_request(:options, METIS_HOST).
    to_return({
    status: 200,
    headers: {
      'Content-Type': "application/json",
    },
    body: route_payload,
  })
  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    .to_return({
      status: 200,
      headers: {
        "Content-Type" => "application/json",
      },
      body: JSON.parse(File.read("spec/fixtures/metis_release_folder_fixture.json")).to_json,
    })

  stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
    .to_return({
      status: 200,
      headers: {
        "Content-Type" => "application/json",
      },
      body: JSON.parse(File.read("spec/fixtures/metis_restrict_folder_fixture.json")).to_json,
    })
end

def stub_list_folder(params = {})
  stub_request(:get, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/#{params[:url_verb] || "list"}\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
      headers: {
        'Content-Type': "application/json",
      },
      body: (params[:response_body] || { files: [], folders: [] }).to_json,
    }).then.to_return({
      status: params[:status_2] || 200,
      headers: {
        'Content-Type': "application/json",
      },
      body: (params[:response_body_2] || { files: [], folders: [] }).to_json,
    })
end

def stub_touch_folder(params = {})
  stub_request(:get, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/touch\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
      headers: {
        'Content-Type': "application/json",
      },
      body: (params[:response_body] || { folders: [] }).to_json,
    })
end

def stub_create_folder(params = {})
  stub_request(:post, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/folder\/create\/#{params[:bucket] || RESTRICT_BUCKET}\//)
    .to_return({
      status: params[:status] || 200,
    })
end

def stub_upload_file(params = {})
  stub_request(:post, /#{METIS_HOST}\/authorize\/upload/)
    .to_return({
      status: params[:status] || 200,
      body: params[:authorize_body] || JSON.generate({}),
    })
  stub_request(:post, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/upload/)
    .to_return({
      status: params[:status] || 200,
      body: params[:upload_body] || JSON.generate({}),
    })
end

def stub_download_file(params = {})
  stub_request(:get, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/download/)
    .to_return({
      status: params[:status] || 200,
      body: params[:file_contents] || "",
    })
end

def stub_bucket_find(params = {})
  stub_request(:post, /#{METIS_HOST}\/#{params[:project] || PROJECT}\/find/)
    .to_return({
      status: params[:status] || 200,
      headers: {
        'Content-Type': 'application/json'
      },
      body: (params[:response_body] || { files: [], folders: [] }).to_json,
    })
end

def stub_magma_models(fixture: "spec/fixtures/magma_test_models.json")
  stub_request(:post, "#{MAGMA_HOST}/retrieve")
    .to_return({ body: File.read(fixture) })
end

def stub_magma_update_json
  stub_request(:post, /#{MAGMA_HOST}\/update$/)
    .to_return({ status: 200, body: {}.to_json, headers: { 'Content-Type': "application/json" } })
end

def stub_redcap_data(stub = nil)
  stub_redcap(
    hash_including(content: "metadata") => File.read("spec/fixtures/redcap_mock_metadata.json"),
    /fields/ => File.read("spec/fixtures/redcap_mock_data_all.json"),
  )
end

def redcap_choices(*choices)
  choices.map.with_index do |c, i|
    "#{i}, #{c}"
  end.join(" | ")
end

def redcap_metadata(template)
  template.map do |form, fields|
    fields.map do |field|
      if field.is_a?(Array)
        field = ["field_name", "field_label", "field_type"].zip(field).to_h.compact
      else
        field = field.map { |k, v| [k.to_s, v] }.to_h
      end

      {
        "field_name": "",
        "form_name": form,
        "section_header": "",
        "field_type": "text",
        "field_label": field["field_name"].split("_").map(&:capitalize).join(" "),
        "select_choices_or_calculations": "",
        "field_note": "",
        "text_validation_type_or_show_slider_number": "",
        "text_validation_min": "",
        "text_validation_max": "",
        "identifier": "",
        "branching_logic": "",
        "required_field": "",
        "custom_alignment": "",
        "question_number": "",
        "matrix_group_name": "",
        "matrix_ranking": "",
        "field_annotation": "",
      }.update(field)
    end
  end.flatten(1).to_json
end

def flat_records(records)
  id_keys = [:record, :redcap_event_name, :redcap_repeat_instrument, :redcap_repeat_instance]

  ret = records.group_by do |r|
    r.values_at(*id_keys)
  end.map do |key_names, eavs|
    id_keys.map(&:to_s).zip(key_names).to_h.compact.merge(
      eavs.map do |eav|
        [eav[:field_name], eav[:value]]
      end.to_h
    ).merge(
      "record_id" => key_names[0],
    )
  end
  ret
end

def redcap_records(base, records)
  records.map do |record|
    base.merge(record)
  end
end

def stub_redcap(stubs)
  stubs.each do |req, rets|
    stub = nil
    [rets].flatten.each do |ret|
      stub = (stub ? stub.then : stub_request(:post, "#{REDCAP_HOST}/api/").with(body: req)).to_return({
        headers: {
          'Content-Type': "application/json",
        },
        body: ret,
      })
    end
  end
end

def stub_redcap_multi_project_records
  stub_redcap(
    /today/ => [
      File.read("spec/fixtures/redcap_mock_data_calendar.json"),
      File.read("spec/fixtures/redcap_mock_data_project_2_calendar.json"),
    ],
    /date_of_birth/ => [
      File.read("spec/fixtures/redcap_mock_data_essential_data.json"),
      File.read("spec/fixtures/redcap_mock_data_project_2_essential_data.json"),
    ],
    /height/ => File.read("spec/fixtures/redcap_mock_data_statistics.json"),
  )
end

def copy_redcap_project(project_name = "test")
  # Make sure the test project is in the right location so the
  # REDCap ETL can find it.
  test_fixture_path = "spec/fixtures/etls/redcap/#{project_name}.rb"
  test_output_path = "#{REDCAP_PROJECT_CONFIG_DIR}/#{project_name}.rb"
  FileUtils.mkdir_p(REDCAP_PROJECT_CONFIG_DIR) unless Dir.exist?(REDCAP_PROJECT_CONFIG_DIR)

  # Make sure we have the newest test project.
  File.delete(test_output_path) if File.file?(test_output_path)
  FileUtils.cp(test_fixture_path, test_output_path)
end

def copy_renaming_project(file_name = "test_renames.json")
  # Make sure the test project is in the right location so the
  # REDCap ETL can find it.
  test_fixture_path = "spec/fixtures/etls/renaming/#{file_name}"
  test_output_path = "#{RENAMING_PROJECT_CONFIG_DIR}/#{file_name}"
  FileUtils.mkdir_p(RENAMING_PROJECT_CONFIG_DIR) unless Dir.exist?(RENAMING_PROJECT_CONFIG_DIR)

  # Make sure we have the newest test config.
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

def stub_ingest_files(file_data = nil)
  file_data ?
    file_data.each do |data|
    create(:ingest_file, **data)
  end : begin
    create(:ingest_file, name: "foo/bar/test1.txt", host: "sftp.example.com", updated_at: "2021-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test2.txt", host: "sftp.example.com", updated_at: "2015-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test3.txt", host: "sftp.example.com", updated_at: "1999-01-01 00:00:00", should_ingest: true)
  end
end

def stub_watch_folders(folder_data = nil)
  folder_data ?
    folder_data.each do |data|
    create(:watch_folder, **data)
  end : begin
    create(:watch_folder, project_name: PROJECT, bucket_name: RELEASE_BUCKET, updated_at: "2021-01-01 00:00:00", folder_path: "path1", watch_type: "link_files", metis_id: 1)
    create(:watch_folder, project_name: PROJECT, bucket_name: RELEASE_BUCKET, updated_at: "2015-01-01 00:00:00", folder_path: "path2", watch_type: "link_files", metis_id: 2)
    create(:watch_folder, project_name: PROJECT, bucket_name: RELEASE_BUCKET, updated_at: "1999-01-01 00:00:00", folder_path: "path1/path1_1", watch_type: "link_files", metis_id: 3)
  end
end

def create_metis_folder(folder_name, folder_path, updated_at: Time.now, id: nil, project_name: PROJECT, bucket_name: RELEASE_BUCKET)
  Etna::Clients::Metis::Folder.new({
    folder_name: folder_name,
    folder_path: folder_path,
    updated_at: updated_at.iso8601,
    id: id,
    project_name: project_name,
    bucket_name: bucket_name
  })
end

def create_metis_file(file_name, file_path, file_hash: SecureRandom.hex, updated_at: Time.now, folder_id: 1)
  Etna::Clients::Metis::File.new({
    file_name: file_name,
    file_path: file_path,
    updated_at: updated_at,
    file_hash: file_hash,
    bucket_name: bucket_name,
    project_name: project_name,
    folder_id: folder_id
  })
end