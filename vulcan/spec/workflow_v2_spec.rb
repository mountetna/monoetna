require_relative '../lib/path'

require 'net/ssh'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:remote_manager) {TestRemoteServerManager.new(Vulcan.instance.ssh_pool)}
  let(:create_workflow_request) {{
    projects: [PROJECT] ,
    repo_url: "/test-utils/available-workflows/snakemake-repo",
    branch: "main",
    workflow_name: "test-workflow",
  }}
  let(:create_workflow_request_all) {{
    projects: ["all"] ,
    repo_url: "/test-utils/available-workflows/snakemake-repo-2", # This does not actually exist
    branch: "test",
    workflow_name: "test-workflow-2",
  }}

  before do
    remove_all_dirs
  end

  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  context 'create workflows' do

    it 'should fail auth if you are not a super user' do
      auth_header(:editor)
      post("/api/v2/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(403)
    end

    it 'should create a workflow object' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)

      # DB object exists
      obj = Vulcan::WorkflowV2.first(name: "test-workflow")
      expect(obj).to_not be_nil
      expect(obj.projects).to eq([PROJECT])
      expect(obj.name).to eq("test-workflow")
      expect(obj.repo_remote_url).to eq("/test-utils/available-workflows/snakemake-repo")
    end

    it 'should inform the user if the workflow exists' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)
      post("/api/v2/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)
      expect(json_body[:msg].include?('exists'))
    end

  end

  context 'list workflows' do

    it 'list workflows available for a project' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      post("/api/v2/workflows/create", create_workflow_request_all)
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workflows].count).to eq(2)
    end

    it 'should return an empty list when no workflows exist' do
      auth_header(:editor)
      get("/api/v2/#{PROJECT}/workflows/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workflows]).to be_empty
    end

  end


  context 'creates workspaces' do

    it 'successfully creates the workspace directory' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        git_version: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}")).to be_truthy
    end

    it 'successfully git clones the workflow' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        git_version: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}/.git")).to be_truthy
    end


    it 'successfully git checkouts the proper tag' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger",
        git_tag: "v1",
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.tag_exists?(obj.path, "v1")).to be_truthy
    end

    it 'successfully uploads the snakemake utils directory' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{obj.path}/snakemake_utils")).to be_truthy
    end

    it 'successfully creates the workspace object' do
      auth_header(:superuser)
      post("/api/v2/workflows/create", create_workflow_request)
      auth_header(:editor)
      request = {
          workflow_id: json_body[:workflow_id],
          workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      expect(last_response.status).to eq(200)
      # DB object exists
      obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(obj).to_not be_nil
      expect(File.basename(obj.path).match?(/\A[a-f0-9]{32}\z/)).to be_truthy
    end

    it 'removes the workspace directory if an error occurs' do
    end

    it 'errors if the requested workflow is not created' do
    end


  end

  context 'list all workspaces' do

    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
    end

    it 'for a user if it exists' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      get("/api/v2/#{PROJECT}/workspace/")
      expect(last_response.status).to eq(200)
    end

    it 'should return an empty list if no workspaces exist' do
    end

    it 'should warn the user if a new version of a workflow exists' do
    end

  end


  context 'list a specific workspace' do
    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)

    end

    it 'for a user if it exists' do
      auth_header(:editor)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
      workspace_id = json_body[:workspace_id]
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
    end

  end

  context 'write files' do

    # TODO: add a file exists endpoint here, we don't want to re-run if the files exist

    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'writes a file to the workspace' do
      file_to_upload = create_temp_file
      auth_header(:editor)
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
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
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

  context 'retrieve dag' do

    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'fetches the dag when no input files need to be generated' do
    end

    it 'generates input files and then fetches the dag' do
      auth_header(:editor)
      workspace_id = Vulcan::Workspace.all[0].id
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/dag")
      expect(json_body[:dag]).to eq(["count", "arithmetic", "checker", "checker_ui", "summary"])
      # Make sure dummy files are removed
    end

  end

  context 'running workflows' do

    # Refer to /spec/fixtures/snakemake-repo/ as the workflow that is being run

    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'can run the first step of a workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # First step in the test workflow is a UI step that writes files to the workspace
      # These are the initial inputs to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
          run: {
            jobs: ["count"],
            params: {
                count_bytes: true,
                count_chars: false
              }
          }
        }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      run_id = json_body[:run_id]
      expect(last_response.status).to eq(200)
      expect(run_id).to_not be_nil

      # Wait until jobs are completed
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end

      # Outputs are created
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to  be_truthy
      # Run objects exist
      obj = Vulcan::Run.first(id: run_id)
      expect(obj).to_not be_nil
      # Correct config file exists
      expect(remote_manager.file_exists?(obj.config_path)).to be_truthy
      config = remote_manager.read_json_file(obj.config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      # Log file exists
      expect(remote_manager.file_exists?(obj.log_path)). to be_truthy
    end

    it 'can run 3 steps of the workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      # Run the first 3 jobs at once
      request = {
        run: {
          jobs: ["count", "arithmetic", "checker"],
          params: {
            count_bytes: true,
            count_chars: false,
            add: 2
          }
        }
      }

      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:run_id]).to_not be_nil
      run_id = json_body[:run_id]

      # Make sure jobs are finished
      check_jobs_status(["count", "arithmetic", "checker"]) do
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
      # Correct config file exists
      expect(remote_manager.file_exists?(obj.config_path)).to be_truthy
      config = remote_manager.read_json_file(obj.config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      expect(config["add"]).to eq("2")
      # Log file exists
      expect(remote_manager.file_exists?(obj.log_path)).to be_truthy
    end


    it 'can run one step and then another' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      # Run the first job
      request = {
        run: {
          jobs: ["count"],
          params: {
            count_bytes: true,
            count_chars: false,
          }
        }
      }

      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end

      # Run the next job
      request = {
        run: {
          jobs: ["arithmetic"],
          params: {
            add: 2
          }
        }
      }
      # Sometimes snakemake still needs a minute to shut-down even though slurm reports the job as complete
      run_workflow_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      end
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]

      # Make sure jobs have finished
      check_jobs_status(["arithmetic"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end

      # Make sure files exist
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/arithmetic.txt")).to be_truthy

      # Make sure two run objects exist
      runs = Vulcan::Run.all
      expect(runs.count).to eq(2)

      # Correct config file exists for run 1
      expect(remote_manager.file_exists?(runs[0].config_path)).to be_truthy
      config = remote_manager.read_json_file(runs[0].config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      expect(config["add"]).to be_nil

      # Correct config file exists for run 2
      expect(remote_manager.file_exists?(runs[1].config_path)).to be_truthy
      config = remote_manager.read_json_file(runs[1].config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      expect(config["add"]).to eq("2")

      # Log file exists
      expect(remote_manager.file_exists?(runs[0].log_path)).to be_truthy
      expect(remote_manager.file_exists?(runs[1].log_path)).to be_truthy
      expect(runs[0].log_path).to_not eq(runs[1].log_path)
    end

    it 'runs the entire workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)

      # Run the first 3 jobs
      request = {
        run: {
          jobs: ["count", "arithmetic", "checker"],
          params: {
            count_bytes: true,
            count_chars: false,
            add: 2
          }
        }
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      run_id = json_body[:run_id]
      check_jobs_status(["count", "arithmetic", "checker"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end

      # Next step involves writing another file to the workspace (checker-ui job)
      request = {
        files: [create_temp_file("ui_check")]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", request, 'CONTENT_TYPE' => 'multipart/form-data')

      # Run the last job, there are no params here
      request = {
        run: {
          jobs: ["summary"],
        }
      }
      run_workflow_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      end

      run_id = json_body[:run_id]

      # Make sure jobs have finished
      check_jobs_status(["summary"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end

      # Make sure files exist
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/count_poem_2.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/check.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/arithmetic.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/ui_check.txt")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/summary.txt")).to be_truthy

      # Make sure two run objects exist
      runs = Vulcan::Run.all
      expect(runs.count).to eq(2)

      # Correct config file exists for run 1
      expect(remote_manager.file_exists?(runs[0].config_path)).to be_truthy
      config = remote_manager.read_json_file(runs[0].config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      expect(config["add"]).to eq("2")

      # Correct config file exists for run 2
      expect(remote_manager.file_exists?(runs[1].config_path)).to be_truthy
      config = remote_manager.read_json_file(runs[1].config_path)
      expect(config["count_bytes"]).to eq("true")
      expect(config["count_chars"]).to eq("false")
      expect(config["add"]).to eq("2")

      # Log file exists
      expect(remote_manager.file_exists?(runs[0].log_path)).to be_truthy
      expect(remote_manager.file_exists?(runs[1].log_path)).to be_truthy
      expect(runs[0].log_path).to_not eq(runs[1].log_path)

    end

    it 'alerts if snakemake is still running' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      request = {
        run: {
          jobs: ["count", "arithmetic", "checker"],
          params: {
            count_bytes: true,
            count_chars: false,
            add: 2
          }
        }
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      expect(last_response.status).to eq(429)
      expect(json_body[:error]).to eq("workflow is still running...")
    end

  end

  context 'status checking' do

    before do
      auth_header(:admin)
      post("/api/v2/repo/clone", create_repo_request)
      post("/api/v2/workflow/publish", publish_workflow_request)
      request = {
        workflow_id: json_body[:workflow_id],
        workspace_name: "running-tiger"
      }
      post("/api/v2/#{PROJECT}/workspace/create", request)
    end

    it 'invokes the first step of a workflow' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # First step in the test workflow requires files to be written to the workspace
      file_names = write_files_to_workspace(workspace.id)
      request = {
        run: {
          jobs: ["count"],
          params: {
            count_bytes: true,
            count_chars: false
          }
        }
      }
      # TODO: add a meta key that can switch profiles
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run", request)
      expect(last_response.status).to eq(200)
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:run_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:count]).to eq("NOT STARTED")
    end
  end

end
