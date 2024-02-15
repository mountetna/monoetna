require_relative '../lib/server/controllers/vulcan_v2_controller'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context 'initialization' do

    it 'lists available workflows' do
      auth_header(:guest)

      # For v1 we can:
      # 1. manually create a directory on c4 with git repos ready to be instantiated
      #    - This should return a bunch of git repos that live on c4 and can be initialized
      #    - Dir structure should be /data2/vulcan/workflows/{project}/
      #    - We can also expose workflows to all projects via: /data2/vulcan/workflows/all/
      get("/api/#{PROJECT}/workflow/list")

    end


    it 'instantiates a new workflow and creates a workspace' do
      auth_header(:guest)

      request = {
        repo: "/data2/vulcan/workflows/#{PROJECT}/test_repo",
        target_dir: "/data2/vulcan/#{PROJECT}",
        workflow_name: "test_workflow"
      }

      post("/api/labors/workflow/init", request)

      expected = {
        workspace: "#{request['target_dir']}/#{request[workflow_name]}/",
        repo: "test-repo/",
      }

      # Test to make sure:
      # - Workspace is created: /data2/vulcan/workflows/labors/test_workflow/hash
      # - Repository is cloned inside the workspace: /data2/vulcan/workflows/labors/test_workflow/hash/test_repo
      # - Record is created in the db
      expect(last_response.status).to eq(200)
    end

    it 'creates a new workspace when the workflow has been changed' do
      # When the git has of the repo has changed
      auth_header(:guest)

      request = {
        repo: "/data2/vulcan/workflows/#{PROJECT}/test_repo",
        target_dir: "/data2/vulcan/#{PROJECT}",
        workflow_name: "test_workflow"
      }

      post("/api/labors/workflow/init", request)

      expected = {
        workspace: "#{request['target_dir']}/#{request[workflow_name]}/",
        repo: "test-repo/",
      }

      # Test to make sure:
      # - New workspace is created: /data2/vulcan/workflows/labors/test_workflow/new_hash
      # - New repository is cloned inside the workspace: /data2/vulcan/workflows/labors/test_workflow/new_hash/test_repo
      # - New Record is created in the db
      expect(last_response.status).to eq(200)

    end
  end

end
