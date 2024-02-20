require_relative '../lib/server/controllers/vulcan_v2_controller'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context 'list pipelines' do

    # For v1 we can:
    # 1. manually create a directory on c4 with git repos ready to be instantiated
    #    - This should return a bunch of git repos that live on c4 and can be initialized
    #    - Dir structure could be /data2/vulcan/pipelines/{project}/
    #    - We can also expose pipelines to all projects via: /data2/vulcan/pipelines/all/

    it 'available for a project' do
      auth_header(:guest)
      get("/api/#{PROJECT}/pipelines/")
      # This should return a list of: :
      # - repo_name
      # - pipeline version
    end

    it 'available for all projects' do
      auth_header(:guest)
      get("/api/all/pipelines/")
    end

  end

  context 'list workspaces' do

    it 'for a pipeline if it exists' do
      pipeline_name = "test_pipeline"
      get("/api/#{PROJECT}/pipeline/#{pipeline_name}")
      # This should return a list of:
      # - workspace hashes
      # - workspace last_updated_at
      # - pipeline version
    end

    it 'should return an empty list if no workspaces exist' do
      pipeline_name = "test_pipeline"
      get("/api/#{PROJECT}/pipeline/#{pipeline_name}")
      # This should return:
      # - workspace hash
      # - workspace last_updated_at
      # - pipeline version
    end

    it 'should warn the user if a new version of a pipeline exists' do
      pipeline_name = "test_pipeline"
      get("/api/#{PROJECT}/pipeline/#{pipeline_name}")
      # This should return:
      # - workspace hash
      # - workspace last_updated_at
      # - pipeline version
      # - A warning that a newer version of a pipeline exists

      # The control flow for this is:
      # - A user has created a workspace and the repo within has a specific version
      # - The repo has been updated on c4
      # - The next time a user does a get request for all pipelines:
      # - There is some logic that determines if their pipeline version is out of date.
      # - The user gets a notification that there is a newer version of the pipeline
      # - They have the option to create a new workspace with that newer version of the pipeline
      # - Note:
      # - When a user first selects a pipeline, it will always be the latest version.
      # - For V1 we should probably not allow them to select versions.
    end

  end

  context 'create workspaces' do

    it 'when the pipeline exists' do
      auth_header(:guest)

      request = {
        repo: "/data2/vulcan/pipelines/#{PROJECT}/test_repo",
        target_dir: "/data2/vulcan/#{PROJECT}",
        pipeline_name: "test_pipeline"
      }

      post("/api/#{PROJECT}/workspace/create", request)

      expected = {
        workspace: "#{request['target_dir']}/#{request[pipeline_name]}/",
        repo: "test-repo/",
      }

      # This should return:
      # - A hash of the workspace that is created
      # This should:
      # - Create a workspace : /data2/vulcan/pipelines/labors/test_pipeline/hash
      # - Cloned repo inside the workspace: /data2/vulcan/pipelines/labors/test_pipeline/hash/test_repo
      # - Record is created in the db
      expect(last_response.status).to eq(200)
    end


  end

  context 'pipeline params' do

    it 'parses the config for pipeline params' do
      # Retrieves all the pipeline params for pipeline
    end

    it 'parses the config for snakemake params' do
      # Retrieves all the snakemake params for the pipeline
    end

  end

  context 'running pipelines' do

    it 'invokes up to n steps of the pipeline' do
      auth_header(:guest)

      request = {
        run_until_step: "3",
        pipeline_params: {

        }
      }
      pipeline_name = "test_pipeline"

      post("/api/#{PROJECT}/#{pipeline_name}/#{workspace_id}/run", request)

      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end

    it 'runs an entire pipeline' do
      # This should return:
      # - A list of output for each job that has run
      # - Whether the job was successful
    end

    it 'reruns a successful job with different pipeline parameters' do
      # This should:
      # - Make sure we mirror the new intermediaries created by snakemake to metis or the Vulcan cache
    end

    it 'reruns a successful job with different snakemake parameters' do
      # this should:
      # - make sure we mirror the new intermediaries created by snakemake to metis or the vulcan cache
    end

    it 'reruns a failed job with different pipeline parameters' do
      # This should:
      # - Make sure we mirror the new intermediaries created by snakemake to metis or the Vulcan cache
    end

    it 'reruns a failed job with different snakemake parameters' do
      # this should:
      # - make sure we mirror the new intermediaries created by snakemake to metis or the vulcan cache
    end

  end

end
