require 'yaml'
require 'logger'
require 'rack/test'
require 'webmock/rspec'
require 'simplecov'
require 'fileutils'
SimpleCov.start
require 'bundler'
Bundler.require(:default, :test)

ENV['METIS_ENV'] = 'test'

require_relative '../lib/metis'
require_relative '../lib/server'
METIS_CONFIG=YAML.load(File.read("config.yml"))
Metis.instance.configure(METIS_CONFIG)
OUTER_APP = Rack::Builder.new do
  use Rack::Static, urls: ['/css', '/js', '/fonts', '/img'], root: 'lib/client'
  use Etna::ParseBody
  use Etna::SymbolizeParams
  use Etna::TestAuth
  use Metis::SetUid
  run Metis::Server.new
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

  config.order = :random

  config.around(:each) do |example|
    DatabaseCleaner.cleaning do
      example.run
    end
  end

  config.before(:suite) do
    stubs.ensure
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
  @stubs ||= Stubs.new
end

def safe_path(path)
  File.join(
    path.split(%r!/!).map do |name|
      Metis::File.safe_file_name(name)
    end
  )
end

class Stubs
  def initialize
    @stubs = []
    @bucket_name = 'files'
  end

  def create_file(project_name, name, contents)
    file_path = project_path(project_name,bucket_path(name))
    stub_file(file_path, contents)
    add_stub(file_path)
  end

  def create_partial(project_name, name, contents, metis_uid)
    partial_path = project_path(
      project_name,
      "uploads/#{Metis::File.safe_file_name("#{metis_uid}-#{name}")}"
    )
    stub_file(partial_path, contents)
    add_stub(partial_path)
  end

  def create_data(project_name, name, contents)
    data_path = project_path(project_name, name)
    stub_file(data_path, contents)
    add_stub(data_path)
  end

  def create_folder(project_name, name)
    folder_path = project_path(project_name, bucket_path(name))
    stub_dir(folder_path)
    add_stub(folder_path)
  end

  def add_file(project_name, name)
    file_path = project_path(project_name, bucket_path(name))
    add_stub(file_path)
  end

  def add_folder(project_name, name)
    folder_path = project_path(project_name, bucket_path(name))
    add_stub(folder_path)
  end

  private

  def bucket_path(name)
    "#{@bucket_name}/#{safe_path(name)}"
  end

  def project_path(project_name, name)
    ::File.expand_path("spec/#{project_name}/#{name}")
  end

  def stub_dir(path)
    FileUtils.mkdir_p(path)
  end

  def add_stub(path)
    return nil if path =~ %r!/files$! || path =~ %r!/\.$!
    stubs.push(path) unless stubs.include?(path)
    return path
  end

  def stub_file(path, contents)
    stub_dir(File.dirname(path))
    File.open(path,"w") do |f|
      f.print contents
    end
  end

  public

  def contents(project)
    [ :athena, :labors ].map do |project|
      [ :uploads, :files ].map do |bucket|
        Dir.glob("spec/#{project}/#{bucket}/*").to_a
      end
    end.flatten
  end

  def ensure
    [ :athena, :labors ].each do |project|
      [ :uploads, :files ].each do |bucket|
        dir = "spec/#{project}/#{bucket}"
        FileUtils.rm_r(dir) if Dir.exists?(dir)
        FileUtils.mkdir_p(dir)
      end
    end
  end

  def clear
    existing_stub_files.each { |stub| File.delete(stub) }
    existing_stub_dirs.each { |stub| FileUtils.rm_r(stub) }
    @stubs = []
  end

  private

  def existing_stub_files
    @stubs.select do |stub|
      File.exists?(stub) && !File.directory?(stub)
    end.sort_by(&:size).reverse
  end

  def existing_stub_dirs
    @stubs.select do |stub|
      File.exists?(stub) && File.directory?(stub)
    end.sort_by(&:size).reverse
  end
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
  },
  viewer: {
    email: 'athena@olympus.org', first: 'Athena', perm: 'v:athena'
  },
  admin: {
    email: 'zeus@olympus.org', first: 'Zeus', perm: 'a:athena'
  }
}
def token_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
end

def hmac_header
  header(*Etna::TestAuth.hmac_header({}))
end

def default_bucket(project_name)
  @default_bucket ||= {}
  @default_bucket[project_name] ||= create( :bucket, project_name: project_name, name: 'files', owner: 'metis', access: 'viewer')
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
