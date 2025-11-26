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
require_relative '../lib/commands'

require 'etna/spec/event_log'

METIS_CONFIG=YAML.load(File.read("config.yml"))
Metis.instance.configure(METIS_CONFIG)

METIS_HOST="metis.#{Metis.instance.config(:token_domain)}"
METIS_URL="https://#{METIS_HOST}"

module Rack::Test::Methods
  def build_rack_mock_session
    Rack::MockSession.new(app, METIS_HOST)
  end
end

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

    # DatabaseCleaner.strategy = :transaction
    DatabaseCleaner.clean_with(:truncation)
  end

  config.order = :random

  config.around(:each) do |example|
    # Unfortunately, DatabaseCleaner + Sequel does not properly handle the auto_savepointing, which means that
    # exceptions handled in rescue blocks do not behave correctly in tests (where as they would be fine outside of
    # tests).  Thus, we are forced to manually handle the transaction wrapping of examples manually to set this option.
    # See: http://sequel.jeremyevans.net/rdoc/files/doc/testing_rdoc.html#label-rspec+-3E-3D+2.8
    #      https://github.com/jeremyevans/sequel/issues/908#issuecomment-61217226
    Metis.instance.db.transaction(:rollback=>:always, :auto_savepoint=>true){ example.run }
  end

  config.before(:suite) do
    stubs.ensure
  end
end

# Helper method to set ledger_tracked_mode_enabled in config
def set_ledger_enabled(value)
  Metis.instance.instance_variable_get(:@config)[:test][:ledger_tracked_mode_enabled] = value
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
  factory :data_block, class: Metis::DataBlock do
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
  post(endpoint.split('/').map { |c| Etna::Route.encode_path_component(c) }.join('/'), hash.to_json, {'CONTENT_TYPE'=> 'application/json'})
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

def glacier_stub(vault_name)
  vault_json = {
      "CreationDate":"2018-08-27T20:56:30.058Z",
      "LastInventoryDate":nil,
      "NumberOfArchives":0,
      "SizeInBytes":0,
      "VaultARN":"arn:aws:glacier:us-east-1:1999:vaults/#{vault_name}",
      "VaultName": vault_name
  }.to_json
  stub_glacier_method(:get, '', body: {Marker: nil,VaultList: []}.to_json)
  stub_glacier_method(:put, "/#{vault_name}", status: 201, body: '{}')
  stub_glacier_method(:get, "/#{vault_name}", body: vault_json)
  stub_glacier_method(:post, "/#{vault_name}", body: vault_json)
  stub_glacier_method(:post, "/#{vault_name}/multipart-uploads",
    status: 201,
    headers: {
      'X-Amz-Multipart-Upload-Id': 'uploadid',
      'Content-Type': 'application/json'
    },
    body: "{}"
  )
  stub_glacier_method(:put, "/#{vault_name}/multipart-uploads/uploadid",
    status: 204
  )
  stub_glacier_method(:post, "/#{vault_name}/multipart-uploads/uploadid",
    status: 201,
    headers: {
      'X-Amz-Archive-Id': 'archive_id',
      'Content-Type': 'application/json',
    }
  )
  stub_glacier_method(:delete, "/#{vault_name}/archives/archive_id",
    status: 204,
    headers: {
      'Content-Type': 'application/json',
    }
  )
end

def stub_glacier_method(method, path, response={})
  stub_request(method, "https://glacier.us-east-1.amazonaws.com/-/vaults#{path}")
    .to_return({
      status: 200,
      headers: {
      'Content-Type' => 'application/json'
      }
    }.merge(response))
end

def stub_metis_download(path, file_data)
  stub_request(:get, /#{path}/)
    .to_return({
      status: 200,
      body: file_data
    })
end

def stub_metis_list_buckets
  return_values = []

  return_values << { status: 502 }
  return_values << {
    status: 200,
    body: {buckets: [{
      bucket_name: 'test',
      project_name: 'athena',
      access: 'viewer'
    }]}.to_json
  }

  stub_request(:get, /list$/)
    .to_return(*return_values)
end

def stub_metis_rm_folder
  return_values = []

  return_values << { status: 502 }
  return_values << {
    status: 200,
    body: {folder_name: ''}.to_json
  }

  stub_request(:delete, /folder\/remove/)
    .to_return(*return_values)
end

def stub_metis_create_folder
  return_values = []

  return_values << { status: 502 }
  return_values << {
    status: 200,
    body: {files: [], folders: []}.to_json
  }

  stub_request(:post, /folder\/create/)
    .to_return(*return_values)
end

def stub_metis_list_bucket(with_error=false)
  stub_request(:get, /list\/bucket/)
    .to_return({
      status: 200,
      body: {files: [], folders: []}.to_json
    })
end

class Stubs
  def initialize
    @stubs = []
  end

  def create_partial(upload, contents, metis_uid)
    partial_path = upload.partial_location
    stub_file(partial_path, contents)
    add_stub(partial_path)
  end

  def create_data(project_name, name, contents)
    data_path = project_path(project_name, name)
    stub_file(data_path, contents)
    add_stub(data_path)
  end

  def create_bucket(project_name, bucket_name)
    bucket_path = project_path(project_name, "buckets/#{bucket_name}")
    stub_dir(bucket_path)
    add_stub(bucket_path)
  end

  def create_folder(project_name, bucket_name, name)
    #folder_path = project_path(project_name, bucket_path(bucket_name, name))
    #stub_dir(folder_path)
    #add_stub(folder_path)
  end

  def create_file(project_name, bucket_name, name, contents, md5_hash=nil)
    hash = md5_hash || Digest::MD5.hexdigest(contents)
    file_path = ::File.expand_path("#{Metis.instance.config(:data_path)}/data_blocks/#{hash[0]}/#{hash[1]}/#{hash}")
    stub_file(file_path, contents)
    add_stub(file_path)
  end

  # have stubs track a file created by someone else
  def add_file(project_name, bucket_name, name)
    file_path = project_path(project_name, bucket_path(bucket_name, name))
    add_stub(file_path)
  end

  def add_folder(project_name, bucket_name, name)
    folder_path = project_path(project_name, bucket_path(bucket_name, name))
    add_stub(folder_path)
  end

  private

  def bucket_path(bucket_name, name)
    "buckets/#{bucket_name}/#{safe_path(name)}"
  end

  def project_path(project_name, name)
    ::File.expand_path("#{Metis.instance.config(:data_path)}/#{project_name}/#{name}")
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
      [ :uploads, :buckets ].each do |bucket|
        dir = "spec/#{project}/#{bucket}"
        FileUtils.rm_r(dir) if Dir.exists?(dir)
        FileUtils.mkdir_p(dir)
      end
    end
  end

  def clear(project_name=nil)
    existing_stub_files.each { |stub| File.delete(stub) }
    existing_stub_dirs.each { |stub| FileUtils.rm_r(stub) }
    FileUtils.rm_r(Dir["#{Metis.instance.config(:data_path)}/#{project_name}/*"]) unless project_name.nil?
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
MAGMA_WISDOM=<<EOT
Although they are
only breath, words
which I command
are immortal
-- (c) Magma
EOT
BACKUP_WISDOM=<<EOT
Although they are
only breath, words
which I command
are immortal
-- copy
EOT
HELMET=<<EOT
  xXx
 xO|Ox
EOT
SHINY_HELMET=<<EOT
 *****
  xXx
 xO|Ox
EOT
LABORS_LIST=<<EOT
The Twelve Labors of Hercules
- The Nemean Lion
- The Lernaean Hydra
- The Ceryneian Hind
- The Erymanthian Boar
- The Augean Stables
- The Stymphalian Birds
- The Cretan Bull
- The Mares of Diomedes
- The Belt of Hippolyta
- The Cattle of Geryon
- The Golden Apples of the Hesperides
- Cerberus
EOT



AUTH_USERS = {
  supereditor: {
    email: 'vesta@olympus.org', name: 'Vesta', perm: 'a:administration'
  },
  editor: {
    email: 'metis@olympus.org', name: 'Metis', perm: 'e:athena,backup'
  },
  restricted_editor: {
    email: 'metis@olympus.org', name: 'Metis', perm: 'E:athena,backup'
  },
  viewer: {
    email: 'athena@olympus.org', name: 'Athena', perm: 'v:athena'
  },
  admin: {
    email: 'zeus@olympus.org', name: 'Zeus', perm: 'a:athena'
  },
  non_user: {
    email: 'nessus@centaurs.org', name: 'Nessus', perm: ''
  },
  guest: {
    email: 'sinon@troy.org', name: 'Sinon', perm: 'g:athena'
  }
}
def token_header(user_type)
  header(*Etna::TestAuth.token_header(AUTH_USERS[user_type]))
end

def hmac_params(params={})
  Etna::TestAuth.hmac_params(params)
end

def hmac_header(params={})
  Etna::TestAuth.hmac_params({
    id: 'metis',
    signature: 'valid'
  }.merge(params)).each do |name, value|
    header( name.to_s, value )
  end
end

def default_bucket(project_name)
  @default_bucket ||= {}
  @default_bucket[project_name] ||= begin
    stubs.create_bucket(project_name, 'files')
    create( :bucket, project_name: project_name, name: 'files', owner: 'metis', access: 'viewer')
  end
end

def create_read_only_bucket(project_name, bucket_name)
  stubs.create_bucket(project_name, bucket_name)
  create( :bucket, project_name: project_name, name: bucket_name, owner: bucket_name, access: 'viewer')
end

def create_upload(project_name, file_name, uid, params={})
  create( :upload,
    {
      project_name: project_name,
      file_name: file_name,
      author: Metis::File.author(Etna::User.new(AUTH_USERS[:editor])),
      metis_uid: uid,
      bucket: params[:bucket] || default_bucket(project_name),
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
      bucket: params[:bucket] || default_bucket(project_name),
      project_name: project_name,
      author: 'metis|Metis',
      folder_name: folder_name,
    }.merge(params)
  )
end

def create_file(project_name, file_name, contents, params={})
  md5_hash = params.delete(:md5_hash) || Digest::MD5.hexdigest(contents)
  
  # Find or create datablock (reuse if same content already exists)
  data_block = Metis::DataBlock.where(md5_hash: md5_hash).first
  unless data_block
    data_block = create(:data_block,
      description: file_name,
      md5_hash: md5_hash,
      size: contents.length,
    )
  end

  create( :file,
    {
      bucket: params[:bucket] || default_bucket(project_name),
      project_name: project_name,
      author: 'metis|Metis',
      file_name: file_name,
      data_block: data_block
    }.merge(params)
  )
end

def upload_file_via_api(project_name, file_name, contents, bucket_name: 'files')
  # Start upload (uses HMAC authentication)
  hmac_header
  json_post(
    "/#{project_name}/upload/#{bucket_name}/#{file_name}",
    file_size: contents.length,
    action: 'start',
    next_blob_size: contents.length,
    next_blob_hash: Digest::MD5.hexdigest(contents)
  )
  expect(last_response.status).to eq(200)
  
  # Upload the file content as a single blob
  blob_file = Tempfile.new('upload_blob')
  blob_file.write(contents)
  blob_file.rewind
  
  hmac_header
  post(
    "/#{project_name}/upload/#{bucket_name}/#{file_name}",
    action: 'blob',
    next_blob_size: 0,
    next_blob_hash: '',
    current_byte_position: 0,
    blob_data: Rack::Test::UploadedFile.new(blob_file.path, 'application/octet-stream')
  )
  expect(last_response.status).to eq(200)
  
  blob_file.close
  blob_file.unlink
  
  # Return the created file (reload to get latest data)
  file = Metis::File.where(file_name: file_name, project_name: project_name).first
  file.reload if file
  
  # Run ChecksumFiles command to compute actual hash from temporary hash
  cmd = Metis::ChecksumFiles.new
  cmd.execute
  file.reload
  
  file
end

def multi_project_backfilled_lifecycle
  wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
  helmet_file = create_file('athena', 'helmet.jpg', HELMET)
  stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
  stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

  labors_file = create_file('labors', 'labors.txt', WISDOM)
  stubs.create_file('labors', 'files', 'labors.txt', WISDOM)

  another_labors_file = create_file('labors', 'labors2.txt', "lol")
  stubs.create_file('labors', 'files', 'labors2.txt', "lol")

  backup_file = create_file('backup', 'backup.txt', WISDOM)
  stubs.create_file('backup', 'files', 'backup.txt', WISDOM)

  wisdom_datablock_id = wisdom_file.data_block_id
  helmet_datablock_id = helmet_file.data_block_id
  labors_datablock_id = labors_file.data_block_id

  # Delete files using the API (simulating real user deletion)
  token_header(:supereditor)
  delete("/athena/file/remove/files/wisdom.txt")
  expect(last_response.status).to eq(200)
  delete("/labors/file/remove/files/labors.txt")
  expect(last_response.status).to eq(200)
  delete("/backup/file/remove/files/backup.txt")
  expect(last_response.status).to eq(200)

  # So now we have:
  # - athena: wisdom.txt (orphaned)
  # - athena: helmet.jpg (not orphaned)
  # - labors: labors2.txt (not orphaned)
  # - labors: labors.txt (orphaned)
  # - backup: backup.txt (orphaned)

  backfill_ledger = Metis::BackfillDataBlockLedger.new
  allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
  backfill_ledger.execute(project_name: 'athena', links: true)
  backfill_ledger.execute(project_name: 'labors', links: true)
  backfill_ledger.execute(project_name: 'backup', links: true)
  allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
  backfill_ledger.execute(orphaned: true)

  {
    wisdom_data_block_id: wisdom_datablock_id,
    helmet_data_block_id: helmet_datablock_id,
    labors_data_block_id: labors_datablock_id,
    another_labors_datablock_id: another_labors_file.data_block_id,
    backup_datablock_id: backup_file.data_block_id
  }
end


def athena_backfilled_lifecycle
  wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
  helmet_file = create_file('athena', 'helmet.jpg', HELMET)
  stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
  stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

  wisdom_datablock_id = wisdom_file.data_block_id
  helmet_datablock_id = helmet_file.data_block_id

  # Delete both files using the API (simulating real user deletion)
  token_header(:editor)
  delete("/athena/file/remove/files/wisdom.txt")
  expect(last_response.status).to eq(200)
  delete("/athena/file/remove/files/helmet.jpg")
  expect(last_response.status).to eq(200)

  backfill_ledger = Metis::BackfillDataBlockLedger.new
  allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
  backfill_ledger.execute(project_name: 'athena', links: true)
  allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
  backfill_ledger.execute(orphaned: true)

  {
    wisdom_data_block_id: wisdom_datablock_id,
    helmet_data_block_id: helmet_datablock_id
  }
end

# From Ruby stdlib tests
# https://fossies.org/linux/ruby/test/readline/test_readline.rb
def with_temp_stdio
  Tempfile.create("test_metis_client_stdin") {|stdin|
    Tempfile.create("test_metis_client_stdout") {|stdout|
      yield stdin, stdout
    }
  }
end

def replace_stdio(stdin_path, stdout_path)
  open(stdin_path, "r"){|stdin|
    open(stdout_path, "w"){|stdout|
      orig_stdin = STDIN.dup
      orig_stdout = STDOUT.dup
      orig_stderr = STDERR.dup
      STDIN.reopen(stdin)
      STDOUT.reopen(stdout)
      STDERR.reopen(stdout)
      begin
        Readline.input = STDIN
        Readline.output = STDOUT
        yield
      ensure
        STDERR.reopen(orig_stderr)
        STDIN.reopen(orig_stdin)
        STDOUT.reopen(orig_stdout)
        orig_stdin.close
        orig_stdout.close
        orig_stderr.close
      end
    }
  }
end

def below_admin_roles
  [:editor, :viewer, :guest]
end

def below_editor_roles
  [:viewer, :guest]
end
