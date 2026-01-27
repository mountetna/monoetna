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
require_relative '../lib/server/controllers/vulcan_v2_controller'
require_relative '../lib/snakemake_command'
require_relative '../lib/snakemake_parser'
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
  superuser: {
    email: 'zeus@twelve-labors.org', name: 'Zeus', perm: 'a:administration', exp: Time.now.to_i + 6000 , flags: 'vulcan'
  },
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

  if ENV['RUN_E2E']=='1'
    config.filter_run e2e: true
  else
    config.filter_run_excluding e2e: true
  end

end

FactoryBot.define do
  factory :workflow, class: Vulcan::WorkflowV2 do
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


def below_admin_roles
  [:editor, :viewer, :guest]
end

def below_editor_roles
  [:viewer, :guest]
end

class TestRemoteServerManager < Vulcan::RemoteManager
  # Auxiliary functions to help with testing but not needed for prod
  def initialize(ssh_pool)
    super(ssh_pool)
  end

  def rmdir(dir)
    # Exclude all the safety checks for testing purposes
    command = Shellwords.join(["rm", "-r", "-f", dir])
    invoke_ssh_command(command)
  end

  def tag_exists?(dir, tag)
    command = "cd #{dir} && git tag -l #{tag}"
    result = invoke_ssh_command(command)
    !result[:stdout].strip.empty?
  end

  def cp_file(src, dest)
    command = "cp #{src} #{dest}"
    invoke_ssh_command(command)
  end

  def create_large_file(remote_file_path, size_mb:)
    Vulcan.instance.logger.info("Creating #{size_mb}MB file at #{remote_file_path}...")
    
    # Create parent directory first
    command = build_command.add('mkdir', '-p', File.dirname(remote_file_path))
    invoke_ssh_command(command.to_s)
    
    Vulcan.instance.logger.info("Generating file content with dd...")
    # Use dd to create a file with random content of specified size
    # bs=1M count=size_mb creates a file of size_mb megabytes
    dd_command = "dd if=/dev/urandom of=#{Shellwords.escape(remote_file_path)} bs=1M count=#{size_mb} 2>/dev/null"
    invoke_ssh_command(dd_command, timeout: 30)
    
    Vulcan.instance.logger.info("Large file created successfully")
    remote_file_path
  end

end



def poem_1_text
  <<~TEXT
    In the realm of the midnight sky,
    Where stars whisper and comets fly,
    A moonlit dance, a celestial show,
    Unfolding secrets we yearn to know.
  TEXT
end

def poem_2_text
  <<~TEXT
    A brook babbles secrets to the stones,
    Tales of ancient earth, of forgotten bones.
    Sunbeams filter through the emerald canopy,
    Painting dappled dreams, a verdant tapestry.
  TEXT
end

def write_files_to_workspace(workspace_id)
  # The first step in the test workflow involves the UI writing files to the workspace
  auth_header(:editor)
  request = {
    files: [{
      filename: "poem.txt",
      content: poem_1_text
    }, {
      filename: "poem_2.txt",
      content: poem_2_text
    }]
  }
  post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request)
  expect(last_response.status).to eq(200)
end

def write_nested_file_to_workspace(workspace_id)
  # The first step in the test workflow involves the UI writing files to the workspace
  auth_header(:editor)
  request = {
    files: [{
      filename: "nested_dir/poem.txt",
      content: poem_1_text
    }]
  }
  post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request)
  expect(last_response.status).to eq(200)
end

def write_image_to_workspace(workspace_id)
  auth_header(:editor)
  request = {
    files: [{
      filename: "image.png",
      content: "image content"
    }]
  }
  post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request)
  expect(last_response.status).to eq(200)
end

def remove_all_dirs
  remote_manager.rmdir(Vulcan::Path.workspace_base_dir)
end


def check_jobs_status(job_names, max_attempts = 5, base_delay = 10)
  attempts = 0

  loop do
    attempts += 1
    yield
    # Check the status of each job in the response
    all_jobs_completed = job_names.all? do |job_name|
      json_body[job_name.to_sym] == "COMPLETED"
    end
    Vulcan.instance.logger.info("job status: #{json_body}")

    # Break the loop if all jobs are completed
    break if all_jobs_completed

    # Break the loop if maximum attempts have been reached
    if attempts >= max_attempts
      raise "Timeout: Maximum attempts reached without all jobs being completed"
    end

    # Sleep with exponential backoff
    sleep_duration = base_delay * (2 ** (attempts - 1))
    sleep(sleep_duration)
  end
end

def run_with_retry(max_attempts = 5, base_delay = 15)
  attempts = 0

  loop do
    attempts += 1
    yield

    Vulcan.instance.logger.info("last response status: #{last_response.status}")
    if last_response.status == 429 || last_response.status == 422 # TODO:revisit - a general retry mechanism?
      if attempts < max_attempts
        sleep_duration = base_delay * (2 ** (attempts - 1))
        sleep(sleep_duration)
      else
        raise "Request failed after #{attempts} attempts due to 429 status."
      end
    else
      break
    end
  end
end

def retry_until(max_attempts: 5, base_delay: 2, &block)
  attempts = 0

  loop do
    attempts += 1
    result = yield
    Vulcan.instance.logger.info("last response: #{json_body}") 
    if result
      break
    end
    
    if attempts >= max_attempts
      raise "Timeout: Condition not met after #{max_attempts} attempts"
    end
    
    # Exponential backoff
    sleep_duration = base_delay * (2 ** (attempts - 1))
    sleep(sleep_duration)
  end
end


def stub_generate_token(project_name, token_name = 'stubbed_token')
  stub_request(:post, /#{Vulcan.instance.config(:janus)[:host]}\/api\/tokens\/generate/)
    .with(
      body: {
        token_type: 'task',
        project_name: project_name,
        read_only: true
      }
    )
    .to_return({
      status: 200,
      body: token_name
    })
end
