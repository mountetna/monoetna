require 'yaml'
require 'logger'
require 'rack/test'
require 'simplecov'
require 'fileutils'
SimpleCov.start
require 'bundler'
Bundler.require(:default, :test)

ENV['METIS_ENV'] = 'test'

require_relative '../lib/metis'
require_relative '../lib/server'

OUTER_APP = Rack::Builder.new do
  use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'lib/client'
  use Etna::ParseBody
  use Etna::SymbolizeParams

  use Etna::TestAuth
  use Metis::SetUid
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

  config.before(:suite) do
    ensure_stub_dirs
  end
end

FactoryBot.define do
  factory :file, class: Metis::File do
    to_create(&:save)
  end
  factory :folder, class: Metis::Folder do
    to_create(&:save)
  end
  factory :bucket, class: Metis::Bucket do
    to_create(&:save)
  end
  factory :upload, class: Metis::Upload do
    to_create(&:save)
  end
end

def fixture(name)
  File.join(File.dirname(__FILE__), "fixtures/#{name}.txt")
end

def json_body
  JSON.parse(last_response.body, symbolize_names: true)
end

def json_post(endpoint, hash)
  post(URI.encode("/#{endpoint}"), hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
end

def stubs
  @stubs ||= []
end

def safe_path(path)
  File.join(
    path.split(%r!/!).map do |name|
      Metis::File.safe_file_name(name)
    end
  )
end

def stub_file(name, contents, project_name = :stub, bucket_name = 'files')
  file_name = "#{bucket_name}/#{safe_path(name)}"
  make_stub(file_name, contents, project_name)
end

def stub_folder(name, project_name = :stub, bucket_name = 'files')
  folder_name = "#{bucket_name}/#{safe_path(name)}"
  make_stub_dir(folder_name, project_name)
end

def stub_data(name, contents, project_name = :stub)
  make_stub(name, contents, project_name)
end

def stub_partial(name, contents, project_name = :stub)
  file_name = "uploads/#{Metis::File.safe_file_name("#{@metis_uid}-#{name}")}"
  make_stub(file_name, contents, project_name)
end

def make_stub_dir(name, project_name)
  return if name == 'files' || name == '.'
  folder_name = "spec/#{project_name}/#{name}"
  FileUtils.mkdir_p(folder_name)
  stubs.push(folder_name) unless stubs.include?(folder_name)
end

def make_stub(name, contents, project_name)
  make_stub_dir(File.dirname(name), project_name)
  file_name = "spec/#{project_name}/#{name}"
  File.open(file_name,"w") do |f|
    f.print contents
  end
  stubs.push(file_name) unless stubs.include?(file_name)
  return File.expand_path(file_name)
end

def ensure_stub_dirs
  [ :athena, :labors ].each do |project|
    [ :uploads, :files ].each do |bucket|
      dir = "spec/#{project}/#{bucket}"
      FileUtils.mkdir_p(dir) unless Dir.exists? dir
    end
  end
end

def clear_stubs
  stub_files = stubs.select do |stub|
    File.exists?(stub) && !File.directory?(stub)
  end.sort_by(&:size).reverse
  stub_dirs = stubs.select do |stub|
    File.exists?(stub) && File.directory?(stub)
  end.sort_by(&:size).reverse

  stub_files.each { |stub| File.delete(stub) }
  stub_dirs.each { |stub| FileUtils.rm_r(stub) }
  @stubs = nil
end

WISDOM=<<EOT
Although they are
only breath, words
which I command
are immortal
EOT
HELMET=<<EOT
  xXx
 xO|Ox
EOT

AUTH_USERS = {
  editor: {
    email: 'metis@olympus.org', first: 'Metis', perm: 'e:athena'
  }
}
def auth_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
end

def hmac_header
  header(*Etna::TestAuth.hmac_header({}))
end

def default_bucket(project_name)
  @default_bucket ||= {}
  @default_bucket[project_name] ||= create( :bucket, project_name: project_name, name: 'files' )
end

def create_upload(project_name, file_name, uid, params={})
  create( :upload,
    {
      project_name: project_name,
      file_name: file_name,
      author: Metis::File.author(Etna::User.new(AUTH_USERS[:editor])),
      metis_uid: uid,
      bucket: default_bucket(project_name),
      current_byte_position: 0,
      file_size: 0,
      next_blob_size: -1,
      next_blob_hash: ''
    }.merge(params)
  )
end

def create_folder(project_name, folder_name, params={})
  create( :folder,
    {
      bucket: default_bucket(project_name),
      project_name: project_name,
      author: 'metis|Metis',
      folder_name: folder_name,
    }.merge(params)
  )
end

def create_file(project_name, file_name, contents, params={})
  create( :file,
    {
      bucket: default_bucket(project_name),
      project_name: project_name,
      author: 'metis|Metis',
      file_name: file_name,
      file_hash: contents ? Digest::MD5.hexdigest(contents) : nil
    }.merge(params)
  )
end
