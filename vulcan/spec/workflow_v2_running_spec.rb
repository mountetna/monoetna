require_relative '../lib/path'

require 'net/ssh'
require 'timeout'

describe VulcanV2Controller do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:remote_manager) { TestRemoteServerManager.new(Vulcan.instance.ssh_pool) }

  before do
    stub_generate_token(PROJECT)
    remove_all_dirs
  end

  def create_workflow(params={})
    create(:workflow, {
      project_name: PROJECT ,
      repo_remote_url: "/test-utils/available-workflows/snakemake-repo",
      name: "test-workflow",
      created_at: Time.now,
      updated_at: Time.now
    }.merge(params))
  end

  def setup_workspace(workflow=nil, params={})
    workflow = create_workflow unless workflow

    project_name = params.delete(:project_name) || PROJECT

    auth_header(:superuser)
    request = {
      workflow_id: workflow.id,
      workspace_name: "running-tiger",
      git_request: "v1"
    }.merge(params)
    post("/api/v2/#{project_name}/workspace/create", request)
    expect(last_response.status).to eq(200)
  end

  context 'running workflows' do

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
end