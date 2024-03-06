require_relative '../lib/server/controllers/vulcan_v2_controller'
require 'net/ssh'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end


  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  # Test that clone repos and make them available to our users
  context 'create workflows' do

    # TODO: this should just be for administrators
    it 'should clone a repo to a directory if it does not exist' do
      auth_header(:guest)
      request = {
        repo: "/test-utils/available-workflows/test-repo",
        workflow_name: "test-workflow",
        branch: "main",
        project_name: PROJECT,
        author: "Jane Doe"
      }
      post("api/v2/workflow/create", request)
      expect(last_response.status).to eq(200)
      msg = "Workflow: #{request[:workflow_name]} successfully cloned and created."
      expect(json_body[:info]).to eql(msg)
      #TODO: assert the directory exists
      expect(Vulcan::WorkflowV2.first(project: request[:project_name], workflow_name: request[:workflow_name])).to be_truthy
    end

    it 'should not clone a repo to a directory if it exists' do
    end



  end

  context 'list workflows' do
    it 'list workflows available for a project' do
      auth_header(:guest)
      request = {
        repo: "/test-utils/available-workflows/test-repo",
        workflow_name: "test-workflow",
        branch: "main",
        project_name: PROJECT,
        author: "Jane Doe"
      }
      post("api/v2/workflow/create", request)
      get("api/v2/#{PROJECT}/workflows/")
      puts last_response
      expect(last_response.status).to eq(200)
      puts json_body
    end

    it 'list workflows available for all projects' do
      auth_header(:guest)
      get("/api/all/workflows/")
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

  context 'create workspaces' do

    it 'when the workflow exists' do
      auth_header(:guest)

      request = {
        repo: "/data2/vulcan/workflows/#{PROJECT}/test-repo",
        target_dir: "/data2/vulcan/#{PROJECT}",
        workflow_name: "test_workflow"
      }

      post("/api/#{PROJECT}/workspace/create", request)

      expected = {
        workspace: "#{request['target_dir']}/#{request[workflow_name]}/",
        repo: "test-repo/",
      }

      # This should return:
      # - A hash of the workspace that is created
      # This should:
      # - Create a workspace : /data2/vulcan/workflows/labors/workflow_name/hash
      # - Cloned repo inside the workspace: /data2/vulcan/workflows/labors/workflow_name/hash/test_repo
      # - Record is created in the db
      # - Parses and saves the workflow params in memory
      expect(last_response.status).to eq(200)
    end

  end

  context 'list workspaces' do

    it 'for a workflow if it exists' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
      # This should return a list of:
      # - workspace hashes
      # - workspace last_updated_at
      # - workflow version
    end

    it 'should return an empty list if no workspaces exist' do
      workflow_name = "test_workflow"
      get("/api/#{PROJECT}/workflow/#{workflow_name}")
      # This should return:
      # - workspace hash
      # - workspace last_updated_at
      # - workflow version
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
