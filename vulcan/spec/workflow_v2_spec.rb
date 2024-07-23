require_relative '../lib/server/controllers/vulcan_v2_controller'
require_relative '../lib/path'

require 'net/ssh'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:remote_manager) {Vulcan::RemoteServerManager.new(Vulcan.instance.ssh_pool)}
  let(:create_repo_request) {{
    project_name: PROJECT,
    repo_url: "/test-utils/available-workflows/test-repo",
    branch: "main",
  }}
  let(:publish_workflow_request) {{
    project_name: PROJECT,
    repo_local_path: "#{Vulcan::Path::WORKFLOW_BASE_DIR}/#{PROJECT}/test-repo",
    workflow_name: "test-workflow",
    branch: "main",
    tag: "v1",
    author: "Jane Doe"
    }
  }

  before do
    remove_all_dirs
  end

  def create_temp_file
    file1 = Tempfile.new(['test_file1', '.txt'])
    file1.write("This is a test file, with content 1")
    file1.rewind
    Rack::Test::UploadedFile.new(file1.path, 'text/plain')
  end

  def remove_all_dirs
    remote_manager.rmdir(Vulcan::Path::WORKFLOW_BASE_DIR, Vulcan::Path::ALLOWED_DIRECTORIES)
    remote_manager.rmdir(Vulcan::Path::WORKSPACE_BASE_DIR, Vulcan::Path::ALLOWED_DIRECTORIES)
    remote_manager.rmdir(Vulcan::Path::TMPDIR, Vulcan::Path::ALLOWED_DIRECTORIES)
  end

  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  context 'clone repo' do

    before do
      remove_all_dirs
    end

    it 'should clone a repo to the project directory' do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      expect(last_response.status).to eq(200)
      # Proper dirs are created
      project_dir = "#{Vulcan::Path::WORKFLOW_BASE_DIR}/#{PROJECT}"
      expect(remote_manager.dir_exists?(project_dir)).to be_truthy
      expect(remote_manager.dir_exists?("#{project_dir}/test-repo")).to be_truthy
    end

    it 'should inform the user if the repo already exists' do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/repo/clone", create_repo_request)
      expect(last_response.status).to eq(200)
      expect(json_body[:Warning].include?('exists'))
    end

    it 'auth should fail if you are not an admin' do
      auth_header(:editor)
      post("/api/v2/#{PROJECT}/", create_repo_request)
      get("/api/v2/#{PROJECT}/repo")
      expect(last_response.status).to eq(403)
    end
  end

  context 'delete repo' do
    it 'auth should fail if you are not a super user' do
      auth_header(:editor)
      delete("/api/v2/#{PROJECT}/test-repo")
      expect(last_response.status).to eq(403)
    end
    it 'auth should delete a repo' do
      auth_header(:superuser)
      post("/api/v2/repo/clone", create_repo_request)
      #delete("/api/v2/#{PROJECT}/test-repo")
      #expect(last_response.status).to eq(200)
    end
  end

  context 'list repos' do
    it 'should list repos in a project directory' do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      get("/api/v2/#{PROJECT}/repo")
      expect(last_response.status).to eq(200)
      expect(json_body[:dirs][0]).to eq("test-repo")
    end

    it 'should list no repos when no repos exist' do
      auth_header(:admin)
      post("/api/v2/repo/create", create_repo_request)
      require 'pry'; binding.pry
    end
  end

  context 'publish workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
    end

    # TODO: this should just be for administrators
    it 'should create a workflow if config is valid' do
      auth_header(:guest)
      post("/api/v2/workflow/publish", publish_workflow_request)
      expect(last_response.status).to eq(200)

      # DB object exists
      obj = Vulcan::WorkflowV2.first(project: PROJECT, workflow_name: "test-workflow")
      expect(obj).to_not be_nil
    end

    it 'should not create a workflow if config is not valid' do
    end


  end

  context 'list workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
    end

    it 'list workflows available for a project' do
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
    end

    it 'list workflows available for all projects' do
    end
  end

  context 'creates workspaces' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
    end

    it 'when the workflow exists' do
      auth_header(:guest)
      get("/api/v2/#{PROJECT}/workflows/")
      request = {
          workflow_id: json_body[:workflows][0][:id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      expect(last_response.status).to eq(200)

      # DB object exists
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(obj).to_not be_nil

      # Proper dirs are created
      workspace_project_dir = "#{Vulcan::Path::WORKSPACE_BASE_DIR}/#{PROJECT}"
      expect(remote_manager.dir_exists?(workspace_project_dir)).to be_truthy
      expect(remote_manager.dir_exists?(obj.path)).to be_truthy

    end

  end

  context 'list workspaces' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end


    it 'for a user if it exists' do
      auth_header(:guest)
      get("/api/v2/#{PROJECT}/workspace/")
      expect(last_response.status).to eq(200)
    end

    it 'should return an empty list if no workspaces exist' do
    end

    it 'should warn the user if a new version of a workflow exists' do
    end

  end


  context 'get workspace' do
    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'for a user if it exists' do
      auth_header(:guest)
      workspace_id = json_body[:workspace_id]
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
    end

  end

  context 'write files' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'writes a file to the workspace' do
      file_to_upload = create_temp_file
      auth_header(:guest)
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        files: [file_to_upload]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request, 'CONTENT_TYPE' => 'multipart/form-data')
      expect(last_response.status).to eq(200)
      # Assert file exists
    end

    it 'writes multiple files' do
    end

  end


  context 'read files' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'requests a file from the workspace' do
      # Write a file to the workspace
      temp_file = create_temp_file
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        files: [temp_file]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request, 'CONTENT_TYPE' => 'multipart/form-data')

      # Read that file
      request = {
        file_names: [temp_file.original_filename]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/read", request)
      expect(last_response.status).to eq(200)
    end

  end

  context 'running workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'makes sure the run_config is valid' do
    end

    it 'invokes 1 step of the workflow' do
      auth_header(:guest)
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        #workflow_run_id: None,
          run: {
            count: {
                poem: "/test-utils/test-input/poem.txt",
                poem_2: "/test-utils/test-input/poem_2.txt",
                count_bytes: true,
                count_chars: false
              }
          }
        }
      # TODO: add a meta key that can switch profiles
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run", request)

      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end

    it 'invokes 1 step of the workflow and checks status' do
      auth_header(:guest)
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        run: {
          count: {
            poem: "/test-utils/test-input/poem.txt",
            poem_2: "/test-utils/test-input/poem_2.txt",
            count_bytes: true,
            count_chars: false
          }
        }
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run", request)
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run/#{json_body[:run_id]}")


      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end


    it 'invokes 2 steps of the workflow' do
      auth_header(:guest)
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        run: {
          count: {
            poem: "/test-utils/test-input/poem.txt",
            poem_2: "/test-utils/test-input/poem_2.txt",
            count_bytes: true,
            count_chars: false
          },
          arithmetic: {
            add: true,
            multiply_by: 10
          }
        }
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run", request)
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run/#{json_body[:run_id]}")
      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end


    it 'runs an entire workflow' do
      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end

    it 'reruns a successful job with different workflow parameters' do
      # This should:
      # - Make sure we mirror the new intermediaries created by snakemake to metis or the Vulcan cache
    end

    it 'reruns a successful job with different snakemake parameters' do
      # this should:
      # - make sure we mirror the new intermediaries created by snakemake to metis or the vulcan cache
    end

    it 'reruns a failed job with different workflow parameters' do
      # This should:
      # - Make sure we mirror the new intermediaries created by snakemake to metis or the Vulcan cache
    end

    it 'reruns a failed job with different snakemake parameters' do
      # this should:
      # - make sure we mirror the new intermediaries created by snakemake to metis or the vulcan cache
    end

  end


end
