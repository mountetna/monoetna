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
      expect(ssh.dir_exists?(obj.path)).to be_truthy

    end

  end

  context 'list workspaces' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
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
      post("/api/v2/repo/create", create_repo_request)
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

  context 'running workflows' do

    before do
      auth_header(:guest)
      post("/api/v2/repo/create", create_repo_request)
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

      workflow_name = "test_workflow"

      post("/api/#{PROJECT}/#{workspace_id}/run", request)

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
