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
    repo: "/test-utils/available-workflows/test-repo",
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
  let(:create_workspace_request) {{
      workflow_name: "test-workflow",
  }}

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
    end
  end

  # Test that clone repos and make them available to our users
  context 'publish workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
    end

    # TODO: this should just be for administrators
    # TODO: should we keep track of all the branches here?
    it 'should clone a repo to a directory if it does not exist' do
      auth_header(:guest)
      post("/api/v2/workflow/publish", publish_workflow_request)
      expect(last_response.status).to eq(200)
      require 'pry'; binding.pry

      # DB object exists
      obj = Vulcan::WorkflowV2.first(project: create_workflow_request[:project_name], workflow_name: create_workflow_request[:workflow_name])
      expect(obj).to_not be_nil
      expect(obj.id).to eq(json_body[:workflow_id])
      expect(obj.workflow_name).to eq(json_body[:workflow_name])

      # Proper dirs are created
      project_dir = "#{WORKFLOW_BASE_DIR}/#{PROJECT}"
      expect(ssh.dir_exists?(project_dir)).to be_truthy
      expect(ssh.dir_exists?(obj.repo_local_path)).to be_truthy

    end

    it 'should not clone a repo to a directory if it exists' do
    end


  end

  context 'list workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/workflow/create", create_workflow_request)
    end

    it 'list workflows available for a project' do
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
    end

    it 'list workflows available for all projects' do
    end
  end

  context 'update workflows' do
    it 'should be able to update a repository' do
    end
  end

  context 'delete workflows' do
    it 'specific permissions' do
    end

  end

  context 'creates workspaces' do

    before do
      auth_header(:guest)
      post("api/v2/workflow/create", create_workflow_request)
    end

    it 'when the workflow exists' do
      # TODO: should we create the workflow by name or id?
      # TODO: should we add a version column or is branch sufficient?
      # How do we want to "version" workflows
      auth_header(:guest)
      require 'pry'; binding.pry
      post("/api/v2/#{PROJECT}/workspace/create", create_workspace_request)
      require 'pry'; binding.pry
      expect(last_response.status).to eq(200)

      # DB object exists
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(obj).to_not be_nil

      # Proper dirs are created
      workspace_project_dir = "#{WORKSPACE_BASE_DIR}/#{PROJECT}"
      expect(ssh.dir_exists?(workspace_project_dir)).to be_truthy
      expect(ssh.dir_exists?(obj.workspace_dir)).to be_truthy

      # Params are saved to to DB

    end

  end

  context 'list workspaces' do

    it 'for a user if it exists' do
      auth_header(:guest)
      post("/api/v2/workflow/create", create_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
      }
      post("/api/v2/#{PROJECT}/workspace/create", create_workspace_request)
      get("/api/v2/#{PROJECT}/workspace/", request)
      expect(last_response.status).to eq(200)
    end

    it 'should return an empty list if no workspaces exist' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
    end

    it 'should warn the user if a new version of a workflow exists' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
      # This should return:
      # - workspace hash
      # - workspace last_updated_at
      # - workflow version
      # - A warning that a newer version of a workflow exists

      # The control flow for this is:
      # - A user has created a workspace and the repo within has a specific version
      # - The repo has been updated on c4
      # - The next time a user does a get request for all workflows:
      # - There is some logic that determines if their workflow version is out of date.
      # - The user gets a notification that there is a newer version of the workflow
      # - They have the option to create a new workspace with that newer version of the workflow
      # - Note:
      # - When a user first selects a workflow, it will always be the latest version.
      # - For V1 we should probably not allow them to select versions.
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
