require_relative '../lib/server/controllers/vulcan_v2_controller'
require 'net/ssh'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:ssh) {Vulcan::SSH.new(Vulcan.instance.ssh)}
  let(:create_repo_request) {{
    project_name: PROJECT,
    repo_url: "/test-utils/available-workflows/test-repo",
    branch: "main",
  }}
  let(:publish_workflow_request) {{
    project_name: PROJECT,
    repo_local_path: "#{WORKFLOW_BASE_DIR}/#{PROJECT}/test-repo",
    workflow_name: "test-workflow",
    branch: "main",
    tag: "v1",
    author: "Jane Doe"
    }
  }

  before do
    remove_all_dirs
  end

  def remove_all_dirs
    ssh.rmdir(WORKFLOW_BASE_DIR, ALLOWED_DIRECTORIES)
    ssh.rmdir(WORKSPACE_BASE_DIR, ALLOWED_DIRECTORIES)
    ssh.rmdir(TMPDIR, ALLOWED_DIRECTORIES)
  end

  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  context 'create repo' do
    it 'should clone a repo to the project directory' do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
      expect(last_response.status).to eq(200)
      require 'pry'; binding.pry

      # Proper dirs are created
      project_dir = "#{WORKFLOW_BASE_DIR}/#{PROJECT}"
      expect(ssh.dir_exists?(project_dir)).to be_truthy
      expect(ssh.dir_exists?("#{project_dir}/test-repo")).to be_truthy
    end
  end

  context 'list repos' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
    end

    it 'should list repos in a project directory' do
      get("/api/v2/#{PROJECT}/repo")
      expect(last_response.status).to eq(200)
      expect(json_body[:dirs][0]).to eq("test-repo")
    end
  end

  context 'publish workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
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
      post("/api/v2/repo/create", create_repo_request)
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
      post("/api/v2/repo/create", create_repo_request)
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
      workspace_project_dir = "#{WORKSPACE_BASE_DIR}/#{PROJECT}"
      expect(ssh.dir_exists?(workspace_project_dir)).to be_truthy
      expect(ssh.dir_exists?(obj.workspace_dir)).to be_truthy

    end

  end

  context 'list workspaces' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      get("/api/v2/#{PROJECT}/workflows/")
      request = {
        workflow_id: json_body[:workflows][0][:id]
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end


    it 'for a user if it exists' do
      auth_header(:guest)
      get("/api/v2/#{PROJECT}/workspace/")
      expect(last_response.status).to eq(200)
    end

    it 'should return an empty list if no workspaces exist' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
    end

    it 'should warn the user if a new version of a workflow exists' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
    end

  end


  context 'workflow params' do
    # Needed for the front-end
    before do
      auth_header(:guest)
      post("/api/v2/workflow/create", create_workflow_request)
      post("/api/v2/#{PROJECT}/workspace/create", create_workspace_request)
    end

    it 'returns a list of workflow params' do
      # Retrieves all the workflow params for workflow
    end

    it 'return a list of config for snakemake params' do
      # Retrieves all the snakemake params for the workflow
    end

  end

  context 'running workflows' do

    it 'invokes up to n steps of the workflow' do
      auth_header(:guest)

      request = {
        run_until_step: "3",
        workflow_params: {

        }
      }
      workflow_name = "test_workflow"

      post("/api/#{PROJECT}/#{workflow_name}/#{workspace_id}/run", request)

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
