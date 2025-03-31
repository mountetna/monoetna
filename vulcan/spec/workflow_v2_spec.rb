require_relative '../lib/path'

require 'net/ssh'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:remote_manager) {TestRemoteServerManager.new(Vulcan.instance.ssh_pool)}
  let(:create_workflow_request) {{
    project_name: PROJECT ,
    repo_url: "/test-utils/available-workflows/snakemake-repo",
    workflow_name: "test-workflow",
  }}

  before do
    remove_all_dirs
  end

  def setup_workspace
    auth_header(:superuser)
    post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
    auth_header(:editor)
    request = {
      workflow_id: json_body[:workflow_id],
      workspace_name: "running-tiger",
      branch: "main",
      git_version: "v1"
    }
    post("/api/v2/#{PROJECT}/workspace/create", request)
    expect(last_response.status).to eq(200)
  end

  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  context 'create workflows' do

    it 'should fail auth if you are not a super user' do
      auth_header(:editor)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(403)
    end

    it 'should create a workflow object' do
      auth_header(:superuser)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)

      # DB object exists
      obj = Vulcan::WorkflowV2.first(name: "test-workflow")
      expect(obj).to_not be_nil
      expect(obj.project_name).to eq(PROJECT)
      expect(obj.name).to eq("test-workflow")
      expect(obj.repo_remote_url).to eq("/test-utils/available-workflows/snakemake-repo")
    end

    it 'should inform the user if the workflow exists' do
      auth_header(:superuser)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)
      expect(json_body[:msg].include?('exists'))
    end

  end

  context 'list workflows' do

    it 'list workflows available for a project' do
      auth_header(:superuser)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workflows].count).to eq(1)
    end

    it 'should return an empty list when no workflows exist' do
      auth_header(:editor)
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workflows]).to be_empty
    end

  end


  context 'creates workspaces', long_running: true do

    before do
      auth_header(:superuser)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
    end

    it 'successfully creates the workspace directory' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        branch: "main",
        git_version: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}")).to be_truthy
    end

    it 'successfully git clones the workflow' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        branch: "main",
        git_version: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}/.git")).to be_truthy
    end


    it 'successfully git checkouts the proper tag' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        branch: "main",
        git_version: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.tag_exists?(obj.path, "v1")).to be_truthy
    end

    it 'successfully uploads the snakemake utils directory' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        branch: "main",
        git_version: "v1"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}/snakemake_utils")).to be_truthy
    end

    it 'successfully creates the workspace object' do
      auth_header(:editor)
      request = {
          workflow_id: json_body[:workflow_id],
          workspace_name: "running-tiger",
          branch: "main",
          git_version: "v1"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      expect(last_response.status).to eq(200)
      # DB object exists
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(obj).to_not be_nil
      expect(obj.dag.to_a).to eq(["count", "arithmetic", "checker", "ui_job_one", "ui_job_two", "summary", "ui_summary", "final"])
      expect(File.basename(obj.path).match?(/\A[a-f0-9]{32}\z/)).to be_truthy
      expect(obj.git_version).to eq("v1")
    end

    it 'successfully sends back vulcan_config and dag' do
      auth_header(:editor)
      request = {
          workflow_id: json_body[:workflow_id],
          workspace_name: "running-tiger",
          branch: "main",
          git_version: "v1"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:vulcan_config]).to_not be_nil
      expect(json_body[:dag]).to_not be_nil
    end

    it 'removes the workspace directory if an error occurs' do
    end

    it 'errors if the requested workflow is not created' do
    end


  end

  context 'list all workspaces' do

    before do
      setup_workspace
    end

    it 'for a user if it exists' do
      get("/api/v2/#{PROJECT}/workspace/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workspaces].count).to eq(1)
    end

    it 'should return an empty list if no workspaces exist' do
    end

    it 'should warn the user if a new version of a workflow exists' do
    end

  end

  context 'saving configs', long_running: true do

    before do
      setup_workspace
    end

    it 'creates a config object and files if it does not exist' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
          params: {
            count_bytes: false,
            count_chars: true
          },
          paramsChanged: [],
          uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:config_id]).to_not be_nil

      # Make sure the config object exists
      config = Vulcan::Config.first(id: json_body[:config_id])
      # Make sure the original config file exists and the default config file
      expect(remote_manager.file_exists?(config.path)).to be_truthy

      # Make sure the config file also has the correct default values as well as the new values
      default_config_content = remote_manager.read_json_file(Vulcan::Path.default_snakemake_config(workspace.path))
      expected_config_content = default_config_content.merge(JSON.parse(request[:params].to_json))
      expect(remote_manager.read_json_file(config.path)).to eq(expected_config_content)
    end

    it 'does not create a new config if the same config already exists' do  
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        paramsChanged: [],
        uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:config_id]).to eq(config_id)
      # Make sure the config file exists
      config = Vulcan::Config.first(id: json_body[:config_id])
      expect(remote_manager.file_exists?(config.path)).to be_truthy
      # Check that there is only one config object
      expect(Vulcan::Config.where(workspace_id: workspace.id).count).to eq(1)
    end

    it 'correctly creates new config objects' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      # Create a config and run a job
      request = {
        params: {
          count_bytes: false,
          count_chars: true 
        },
        paramsChanged: [],
        uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_1 = Vulcan::Config.first(id: json_body[:config_id])
      request_2 = {
        params: {
          count_bytes: false,
          count_chars: false 
        },
        paramsChanged: ["count_chars"],
        uiFilesSent: []
      }

      # NOTE: when we set
      # params: {
      #   count_bytes: true,
      #   count_chars: false 
      # }
      # after the first request snakemake does not detect that the config has changed.
      # Not entirely sure why, but this just happens with our example workflow.


      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request_2)
      config_2 = Vulcan::Config.first(id: json_body[:config_id])

      # Make sure the config object exists
      expect(Vulcan::Config.where(workspace_id: workspace.id).count).to eq(2)

      # Make sure the the 2 config files exists
      expect(remote_manager.file_exists?(config_1.path)).to be_truthy
      expect(remote_manager.file_exists?(config_2.path)).to be_truthy

      config_content_1 = remote_manager.read_json_file(config_1.path)
      config_content_2 = remote_manager.read_json_file(config_2.path)

      # Check that the 2 config files are different
      expect(config_content_1).to_not eq(config_content_2)

      # Confirm that the config files have the correct values
      request[:params].each do |key, expected_value|
        expect(config_content_1[key.to_s]).to eq(expected_value)
      end

      request_2[:params].each do |key, expected_value|
        expect(config_content_2[key.to_s]).to eq(expected_value)
      end

    end

    it 'correctly returns scheduled and downstream jobs' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
          params: {
            count_bytes: false,
            count_chars: true,
            add: 2,
            add_and_multiply_by: 2
          },
          paramsChanged: [],
          uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(json_body[:scheduled]).to match_array(["arithmetic", "checker", "count"])
      expect(json_body[:downstream]).to match_array(["ui_job_one", "ui_job_two", "summary", "ui_summary", "final"])
    end
  end

  context 'list a specific workspace', long_running: true do
    before do
      setup_workspace
    end

    it 'returns a config if no run exists' do
      workspace_id = json_body[:workspace_id]
      workspace = Vulcan::Workspace.first(id: workspace_id)
      write_files_to_workspace(workspace_id)
      # Create a config and run a job
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        paramsChanged: [],
        uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/config", request)
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
      # Make sure the config file also has the correct default values as well as the new values
      default_config_content = remote_manager.read_json_file(Vulcan::Path.default_snakemake_config(workspace.path))
      expected_config_content = default_config_content.merge(JSON.parse(request[:params].to_json))
      expect(json_body[:last_config]).to eq(expected_config_content.transform_keys(&:to_sym))
      expect(json_body[:last_config_id]).to_not be_nil
      expect(json_body[:last_job_status]).to be_nil
      expect(json_body[:dag]).to_not be_nil
      expect(json_body[:vulcan_config]).to_not be_nil
      expect(json_body[:last_run_id]).to be_nil
    end

    it 'returns the last run and last config' do
      workspace_id = json_body[:workspace_id]
      write_files_to_workspace(workspace_id)
      # Create a config and run a job
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        paramsChanged: [],
        uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/config", request)
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run/#{run_id}")
      end

      # Create another config and run a job
      request_2 = {
        params: {
          count_bytes: false,
          count_chars: false 
        },
        paramsChanged: ["count_chars"],
        uiFilesSent: []
      }
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/config", request_2)
      end
      config_id = json_body[:config_id]
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/run/#{config_id}")
      end

      # Check that the last config is the second config
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      workspace = Vulcan::Workspace.first(id: workspace_id)
      default_config_content = remote_manager.read_json_file(Vulcan::Path.default_snakemake_config(workspace.path))
      expected_config_content = default_config_content.merge(JSON.parse(request_2[:params].to_json))
      expect(json_body[:last_config]).to eq(expected_config_content.transform_keys(&:to_sym))
      expect(json_body[:last_config_id]).to_not be_nil
      expect(json_body[:last_job_status]).to_not be_nil
      expect(json_body[:last_run_id]).to_not be_nil
      expect(json_body[:dag]).to_not be_nil
      expect(json_body[:vulcan_config]).to_not be_nil
    end

    it 'returns no configs if they do not exist' do
      workspace_id = json_body[:workspace_id]
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:last_config]).to be_nil
      expect(json_body[:last_job_status]).to be_nil
    end


  end

  context 'write files' do

    # TODO: add a file exists endpoint here, we don't want to re-run if the files exist

    before do
      setup_workspace
    end

    it 'writes a file to the workspace' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      file_name = "test_file_name"
      content = "This is a test file, with content 1."
      request = {
        files: [{
          filename: file_name,
          content: content
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", request)
      expect(last_response.status).to eq(200)
      expect(remote_manager.file_exists?("#{workspace.path}/output/#{file_name}")).to be_truthy
      expect(remote_manager.read_file_to_memory("#{workspace.path}/output/#{file_name}").chomp).to eq(content)
    end

    it 'writes multiple files' do
    end

  end


  context 'read files' do

    before do
      setup_workspace
    end

    it 'requests a file from the workspace' do
      # Write a file to the workspace
      workspace_id = Vulcan::Workspace.all[0].id
      request = {
        files: [{
          filename: "test_file_name",
          content: "This is a test file, with content 1"
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request)
      expect(last_response.status).to eq(200)

      # Read that file
      request = {
        file_names: ["test_file_name"]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/read", request)
      expect(last_response.status).to eq(200)
    end

  end

  context 'get files' do

    before do
      setup_workspace
    end

    it 'gets a list of files in the workspace' do
      workspace_id = Vulcan::Workspace.all[0].id
      write_files_to_workspace(workspace_id)
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/")
      expect(last_response.status).to eq(200)
      expect(json_body[:files]).to eq(["poem.txt", "poem_2.txt"])
    end

  end

  context 'update workspace' do
    before do
      setup_workspace
    end

    it 'updates the workspace name and tags' do
      workspace_id = Vulcan::Workspace.all[0].id
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/update", {name: "new_name", tags: ["new_tag", "new_tag_2"]})
      expect(last_response.status).to eq(200)
      workspace = Vulcan::Workspace.first(id: workspace_id)
      expect(workspace.name).to eq("new_name")
      expect(workspace.tags).to eq(["new_tag", "new_tag_2"])
    end

  end

  
  context 'running workflows', long_running: true do

    # Refer to /spec/fixtures/snakemake-repo/ as the workflow that is being run

    before do
      setup_workspace
    end

    def write_file_to_workspace(original_request, workspace, file_name, content)
      file_request = {
        files: [{
          filename: file_name,
          content: content
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", file_request)
      expect(last_response.status).to eq(200)
      config_request = original_request.merge(uiFilesSent: ["output/#{file_name}"])
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", config_request)
      end
    end

    def run_entire_test_workflow(workspace)
      auth_header(:editor)
      write_files_to_workspace(workspace.id)

      request_first_jobs = {
        params: {
          count_bytes: true,
          count_chars: false,
          add: 2,
          add_and_multiply_by: 4
        },
        uiFilesSent: [],
        paramsChanged: []
      }

      # Run the first 3 jobs
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request_first_jobs)
      config_id = json_body[:config_id]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      run_id = json_body[:run_id]
      check_jobs_status(["count", "arithmetic", "checker"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)

      # Next step involves UI steps - so we just simulate this by writing files to the workspace
      write_file_to_workspace(request_first_jobs, workspace, "ui_job_one.txt", "This is a test file, with content 1")
      write_file_to_workspace(request_first_jobs, workspace, "ui_job_two.txt", "This is a test file, with content 2")
      config_id = json_body[:config_id]
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      end
      expect(last_response.status).to eq(200)

      # Make sure jobs have finished
      run_id = json_body[:run_id]
      check_jobs_status(["summary"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)

      write_file_to_workspace(request_first_jobs, workspace, "ui_summary.txt", "This is a test file, with content 3")
      # Run the final job
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      end
      run_id = json_body[:run_id]
      check_jobs_status(["final"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)
    end

    it 'alerts if snakemake is still running' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      request = {
          params: {
            count_bytes: true,
            count_chars: false,
            add: 2,
            add_and_multiply_by: 4
          },
          uiFilesSent: [],
          paramsChanged: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      expect(last_response.status).to eq(429)
      expect(json_body[:error]).to eq("workflow is still running...")
    end


    it 'can run the first step of a workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # First step in the test workflow is a UI step that writes files to the workspace
      # These are the initial inputs to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        uiFilesSent: [],
        paramsChanged: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(last_response.status).to eq(200)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      run_id = json_body[:run_id]
      expect(last_response.status).to eq(200)
      expect(run_id).to_not be_nil
      # Wait until jobs are completed
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      # Outputs are created
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      # Run objects exist
      obj = Vulcan::Run.first(id: run_id)
      expect(obj).to_not be_nil
      # Correct config file exists
      expect(remote_manager.file_exists?(obj.log_path)). to be_truthy
    end

    it 'can run 3 steps of the workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      # Run the first 3 jobs at once
      request = {
        params: {
            count_bytes: false,
            count_chars: true,
            add: 2,
            add_and_multiply_by: 2
        },
        uiFilesSent: [],
        paramsChanged: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:run_id]).to_not be_nil
      run_id = json_body[:run_id]
      # Make sure jobs are finished
      check_jobs_status(["count", "arithmetic", "checker"], 5) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      # Outputs are created
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/arithmetic.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/check.txt")).to be_truthy
      # Run objects exist
      obj = Vulcan::Run.first(id: run_id)
      expect(obj).to_not be_nil
      # Log file exists
      expect(remote_manager.file_exists?(obj.log_path)).to be_truthy
    end

    it 'can run one step and then another' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      # Run the first job
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        uiFilesSent: [],
        paramsChanged: []
      }

      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(last_response.status).to eq(200)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)

      # Run the next job
      request = {
        params: {
          count_bytes: false,
          count_chars: true,
          add: 2,
          add_and_multiply_by: 4
        },
        uiFilesSent: [],
        paramsChanged: []
      }
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      end
      expect(last_response.status).to eq(200)
      config_id = json_body[:config_id]
      # Sometimes snakemake still needs a minute to shut-down even though slurm reports the job as complete
      run_with_retry(5) do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      end
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]
      # Make sure jobs have finished
      check_jobs_status(["arithmetic", "checker"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)

      # Make sure files exist
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/arithmetic.txt")).to be_truthy

      # Make sure two run objects exist
      runs = Vulcan::Run.all
      expect(runs.count).to eq(2)

      # Log file exists
      expect(remote_manager.file_exists?(runs[0].log_path)).to be_truthy
      expect(remote_manager.file_exists?(runs[1].log_path)).to be_truthy
      expect(runs[0].log_path).to_not eq(runs[1].log_path)
    end

    it 'runs the entire workflow' do
      workspace = Vulcan::Workspace.all[0]
      run_entire_test_workflow(workspace)
      expect(last_response.status).to eq(200)

      # Make sure files exist
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/check.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/arithmetic.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/summary.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_summary.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/final.txt")).to be_truthy

      # Make sure run objects exist
      runs = Vulcan::Run.all
      expect(runs.count).to eq(3)

      # Log file exists
      expect(remote_manager.file_exists?(runs[0].log_path)).to be_truthy
      expect(remote_manager.file_exists?(runs[1].log_path)).to be_truthy
      expect(remote_manager.file_exists?(runs[2].log_path)).to be_truthy
      expect(runs[0].log_path).to_not eq(runs[1].log_path)
      expect(runs[1].log_path).to_not eq(runs[2].log_path)
      expect(runs[0].log_path).to_not eq(runs[2].log_path)

    end

    it 'halts execution of a UI step after a param has changed' do
      # Run entire workflow
      workspace = Vulcan::Workspace.all[0]
      run_entire_test_workflow(workspace)
      # Now we go back and change the param of earlier job
      request = {
        params: {
          count_bytes: true,
          count_chars: false,
          add: 4, # Changed from 2
          add_and_multiply_by: 4
        },
        uiFilesSent: [],
        paramsChanged: ["add"]
      }
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      end
      # Make sure that we removed all ui files
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_job_one.txt")).to be_falsey
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_job_two.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_summary.txt")).to be_falsey
    end

    it 'halts execution of a UI step after a UI step has changed' do
      # Run entire workflow
      workspace = Vulcan::Workspace.all[0]
      run_entire_test_workflow(workspace)
      request_first_jobs = {
        params: {
          count_bytes: true,
          count_chars: false,
          add: 2,
          add_and_multiply_by: 4
        },
        uiFilesSent: [],
        paramsChanged: []
      }
      write_file_to_workspace(request_first_jobs, workspace, "ui_job_one.txt", "This is a test file, with content haha")
      # Next step involves UI steps - so we just simulate this by writing files to the workspace
      # Change a file
      ui_file_request = {
        params: {
          count_bytes: true,
          count_chars: false,
          add: 2,
          add_and_multiply_by: 4
        },
        uiFilesSent: ["output/ui_job_one.txt"], # Do we want to include the output/ prefix?
        paramsChanged: []
      }
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", ui_file_request)
      end
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_job_one.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_job_two.txt")).to be_truthy # This should not be removed
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_summary.txt")).to be_falsey
    end
  end

  context 'status checking', long_running: true do

    before do
      setup_workspace
    end

    it 'checks the status of the first step of a workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      request = {
          params: {
            count_bytes: false,
            count_chars: true 
          },
          uiFilesSent: [],
          paramsChanged: []
      }
      # TODO: add a meta key that can switch profiles
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      expect(last_response.status).to eq(200)
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:run_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:count]).to eq("NOT STARTED")
    end
  end

  context 'is running' do

    before do
      setup_workspace
    end

    it 'returns false if snakemake is is not running' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/running")
      expect(last_response.status).to eq(200)
      expect(json_body[:running]).to be_falsey
    end

    it 'returns true if snakemake is running' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      request = {
          params: {
            count_bytes: false,
            count_chars: true 
          },
          paramsChanged: [],
          uiFilesSent: []
      }
      # TODO: add a meta key that can switch profiles
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/running")
      expect(last_response.status).to eq(200)
      expect(json_body[:running]).to be_truthy
    end
  end

  context 'read image' do

    before do
      setup_workspace
    end

    it 'returns the image' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_image_to_workspace(workspace.id)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/image/read", {file_name: "image.png"})
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq("image content\n")
    end

    it 'raises an error if the file does not exist' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/image/read", {file_name: "does_not_exist.png"})
      expect(last_response.status).to eq(422)
    end
    
  end

end
