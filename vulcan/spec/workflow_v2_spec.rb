require_relative '../lib/path'

require 'net/ssh'

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
      branch: "main",
      git_request: "v1"
    }.merge(params)
    post("/api/v2/#{project_name}/workspace/create", request)
    expect(last_response.status).to eq(200)
  end

  context 'ssh' do
    it 'should warn the user and start the API if we cannot establish a ssh connection' do
    end
  end

  context 'create workflows' do
    let(:create_workflow_request) {{
      project_name: PROJECT ,
      repo_url: "/test-utils/available-workflows/snakemake-repo",
      workflow_name: "test-workflow",
    }}

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
      create_workflow

      auth_header(:superuser)
      post("/api/v2/#{PROJECT}/workflows/create", create_workflow_request)
      expect(last_response.status).to eq(200)
      expect(json_body[:msg].include?('exists'))
    end

  end

  context 'list workflows' do

    it 'list workflows available for a project' do
      create_workflow

      auth_header(:viewer)
      get("/api/v2/#{PROJECT}/workflows")

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


  context 'creates workspaces' do
    let(:workflow) { create_workflow }
    let(:workspace_request) {
      {
        workflow_id: workflow.id,
        workspace_name: "running-tiger",
        branch: "main",
        git_request: "v1",
      }
    }

    def post_workspace_create(params)
      post("/api/v2/#{PROJECT}/workspace/create", params)
    end

    it 'successfully creates the workspace directory' do
      auth_header(:editor)
      post_workspace_create(workspace_request)

      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{workspace.path}")).to be_truthy
    end

    it 'successfully git clones the workflow' do
      auth_header(:editor)
      post_workspace_create(workspace_request)

      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.dir_exists?("#{workspace.path}/.git")).to be_truthy
    end

    it 'successfully creates a keep file in the output directory' do
      auth_header(:editor)
      post_workspace_create(workspace_request)

      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.file_exists?("#{workspace.path}/output/.keep")).to be_truthy
    end

    it 'successfully git checkouts the proper tag' do
      auth_header(:editor)
      post_workspace_create(workspace_request)
      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.tag_exists?(workspace.path, "v1")).to be_truthy
    end

    it 'successfully creates the dl_config.yaml file' do
      auth_header(:editor)
      post_workspace_create(workspace_request)
      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(remote_manager.file_exists?("#{workspace.path}/dl_config.yaml")).to be_truthy
      expect(remote_manager.read_yaml_file("#{workspace.path}/dl_config.yaml")).to eq({
        "project_name" => "#{PROJECT}",
        "token" => "stubbed_token",
        "magma_url" => "#{Vulcan.instance.config(:magma)[:host]}",
      })
    end

    # For now we are just using the default profile

    # it 'successfully uploads the profiles directory' do
    #   auth_header(:editor)
    #   request = {
    #     workflow_id: json_body[:workflow_id],
    #     workspace_name: "running-tiger",
    #     branch: "main",
    #     git_request: "v1"
    #   }
    #   post("/api/v2/#{PROJECT}/workspace/create", request)
    #   obj = Vulcan::Workspace.first(id: json_body[:workspace_id])
    #   expect(remote_manager.dir_exists?("#{obj.path}/profiles/generic")).to be_truthy
    # end

    it 'successfully creates the workspace object' do
      auth_header(:editor)
      post_workspace_create(workspace_request)
      expect(last_response.status).to eq(200)
      workspace = Vulcan::Workspace.first(id: json_body[:workspace_id])
      expect(workspace).to_not be_nil
      expect(workspace.dag).to eq({
        "final"=>["ui_summary"],
        "ui_summary"=>["ui_job_one", "ui_job_two", "summary"],
        "ui_job_one"=>["checker"],
        "checker"=>["arithmetic"],
        "arithmetic"=>["count"],
        "count"=>[],
        "ui_job_two"=>[],
        "summary"=>["count", "arithmetic", "checker", "ui_job_one", "ui_job_two"]
      })
      expect(File.basename(workspace.path).match?(/\A[a-f0-9]{32}\z/)).to be_truthy
      expect(workspace.git_ref).to eq("v1")
      expect(workspace.git_sha).to_not eq(nil)
    end

    it 'successfully sends back vulcan_config and dag, file_dag, and dag_flattened' do
      auth_header(:editor)
      post_workspace_create(workspace_request)
      expect(last_response.status).to eq(200)
      expect(json_body[:vulcan_config]).to_not be_nil
      expect(json_body[:dag]).to_not be_nil
      expect(json_body[:dag_flattened]).to_not be_nil
      expect(json_body[:file_dag]).to_not be_nil
    end
  end

  context 'list project workspaces' do
    it 'for a user if it exists' do
      setup_workspace
      auth_header(:editor)
      get("/api/v2/#{PROJECT}/workspace")
      expect(last_response.status).to eq(200)
      expect(json_body[:workspaces].count).to eq(1)
    end

    it 'does not list other projects' do
      setup_workspace

      stub_generate_token('athena')
      workflow2 = create_workflow(name: "other-workflow", project_name: "athena")
      setup_workspace(workflow2, project_name: "athena", workspace_name: "walking-owl")

      auth_header(:editor)
      get("/api/v2/#{PROJECT}/workspace")
      expect(last_response.status).to eq(200)
      expect(json_body[:workspaces].length).to eq(1)
      expect(json_body[:workspaces]).to all(include(:workspace_path => a_string_starting_with("#{Vulcan.instance.config(:base_dir)}/workspace/#{PROJECT}")))
    end

    it 'should return an empty list if no workspaces exist' do
      auth_header(:editor)
      get("/api/v2/#{PROJECT}/workspace/")
      expect(last_response.status).to eq(200)
      expect(json_body[:workspaces]).to be_empty
    end
  end

  context 'saving configs' do

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

    it 'creates a new config if the same config already exists' do  
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
      # Make sure the config file exists
      config = Vulcan::Config.all
      expect(config.count).to eq(2)
      config_1 = config.first
      config_2 = config.last
      expect(remote_manager.file_exists?(config_1.path)).to be_truthy
      expect(remote_manager.file_exists?(config_2.path)).to be_truthy
      expect(config_1.path).to_not eq(config_2.path)
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

      # Confirm that the params are returned in the response 
      expect(json_body[:params]).to eq(request_2[:params])
    end

    it 'correctly returns the scheduled and downstream files (first job)' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first snakemake job
      request = {
          params: {
            count_bytes: false,
            count_chars: true,
          },
          paramsChanged: [],
          uiFilesSent: []
      }

      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(json_body[:jobs][:planned]).to match_array(["count"])
      expect(json_body[:jobs][:unscheduled]).to match_array(["arithmetic", "checker", "ui_job_one", "ui_job_two", "summary", "ui_summary", "final"])
      expect(json_body[:jobs][:completed]).to match_array([])

      expect(json_body[:files][:planned]).to match_array(["output/count_poem.txt", "output/count_poem_2.txt"])
      expect(json_body[:files][:unscheduled]).to match_array(["output/arithmetic.txt", "output/check.txt", "output/ui_job_one.txt", "output/ui_job_two.txt", "output/summary.txt", "output/ui_summary.txt", "output/final.txt"])
      expect(json_body[:files][:completed]).to match_array([])
       
      config = Vulcan::Config.first(id: json_body[:config_id])
      expect(config.input_files).to eq(["output/poem.txt", "output/poem_2.txt", "resources/number_to_add.txt"])
      expect(config.input_params.to_h.transform_values(&:to_s)).to eq(request[:params].transform_keys(&:to_s).transform_values(&:to_s))
    end


    it 'correctly returns scheduled and downstream files and jobs (three jobs)' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first 3 snakemake jobs
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

      # TODO: figure maybe change algorithm to include ui_job_two and output/ui_job_two.txt
      # Both are root nodes but should probably be included
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      expect(json_body[:jobs][:completed]).to match_array([])
      expect(json_body[:jobs][:planned]).to match_array(["count", "arithmetic", "checker"])
      expect(json_body[:jobs][:unscheduled]).to match_array(["ui_job_one", "ui_job_two", "summary", "ui_summary", "final"])

      expect(json_body[:files][:completed]).to match_array([])
      expect(json_body[:files][:planned]).to match_array(["output/count_poem.txt", "output/count_poem_2.txt", "output/arithmetic.txt", "output/check.txt"])
      expect(json_body[:files][:unscheduled]).to match_array(["output/ui_job_one.txt", "output/ui_job_two.txt", "output/summary.txt", "output/ui_summary.txt", "output/final.txt"])

      config = Vulcan::Config.first(id: json_body[:config_id])
      expect(config.input_files).to eq(["output/poem.txt", "output/poem_2.txt", "resources/number_to_add.txt"])
      expect(config.input_params.to_h.transform_values(&:to_s)).to eq(request[:params].transform_keys(&:to_s).transform_values(&:to_s))
    end

    it 'it correctly returns scheduled and downstream files and jobs after a config has been run and then changed' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      # We need to write some initial input files to the workspace.
      write_files_to_workspace(workspace.id)
      # Next we run the first 3 snakemake jobs
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
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      run_id = json_body[:run_id]
      check_jobs_status(["count", "arithmetic", "checker"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}") 
      end 
      expect(last_response.status).to eq(200)

      # Now we change the config and run the workflow again
      request_2 = {
        params: {
          count_bytes: false,
          count_chars: true,
          add: 4, # Change the add param
          add_and_multiply_by: 2
        },
        paramsChanged: ["add"],
        uiFilesSent: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request_2)
      config_request = json_body
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:config_id]}")
      run_id = json_body[:run_id]
      check_jobs_status(["arithmetic", "checker"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}") 
      end 
      expect(last_response.status).to eq(200)
      expect(config_request[:jobs][:completed]).to match_array(["count"])
      expect(config_request[:jobs][:planned]).to match_array(["arithmetic", "checker"])
      expect(config_request[:jobs][:unscheduled]).to match_array(["ui_job_one", "ui_job_two", "summary", "ui_summary", "final"])

      expect(config_request[:files][:completed]).to match_array(["output/count_poem.txt", "output/count_poem_2.txt"])
      expect(config_request[:files][:planned]).to match_array(["output/arithmetic.txt", "output/check.txt"])
      expect(config_request[:files][:unscheduled]).to match_array(["output/ui_job_one.txt", "output/ui_job_two.txt", "output/summary.txt", "output/ui_summary.txt", "output/final.txt"])

    end

  end

  context 'list a specific workspace' do
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
      expect(json_body[:dag_flattened]).to_not be_nil
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

    it 'updates the token in the dl_config.yaml' do
      workspace_id = json_body[:workspace_id]
      stub_generate_token(PROJECT, "another_token")
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
      expect(remote_manager.read_yaml_file(Vulcan::Path.dl_config(json_body[:workspace_path]))["token"]).to eq("another_token")
    end

  end

  context 'write files' do

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
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      file_name = "test_file_name"
      content = "This is a test file, with content 1."
      file_name_2 = "test_file_name_2"
      content_2 = "This is a test file, with content 2."
      request = {
        files: [{
          filename: file_name,
          content: content
        }, {
          filename: file_name_2,
          content: content_2
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", request)
      expect(last_response.status).to eq(200)
      expect(remote_manager.file_exists?("#{workspace.path}/output/#{file_name}")).to be_truthy
      expect(remote_manager.file_exists?("#{workspace.path}/output/#{file_name_2}")).to be_truthy
    end

    it 'handles a json string correctly' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      file_name = "test_file_name"
      json_content = '{"contents":[{"col":"gene1","def":["exactly",0,"below",0.02],"logic":null}]}'
      request = {
        files: [{
          filename: file_name,
          content: json_content
        }]
      }
      header 'Content-Type', 'application/json'
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", request.to_json)
      expect(last_response.status).to eq(200)
      expect(remote_manager.file_exists?("#{workspace.path}/output/#{file_name}")).to be_truthy
      expect(remote_manager.read_file_to_memory("#{workspace.path}/output/#{file_name}").chomp).to eq(json_content)
    end

  end


  context 'read files' do

    before do
      setup_workspace
    end

    it 'requests a UTF-8 file from the workspace' do
      workspace_id = Vulcan::Workspace.all[0].id
      utf8_content = "This is a test file, with UTF-8 content: こïv"
      request = {
        files: [{
          filename: "utf8_file.txt",
          content: utf8_content
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/write", request)
      expect(last_response.status).to eq(200)

      request = {
        file_names: ["utf8_file.txt"]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace_id}/file/read", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:files][0][:content]).to eq(utf8_content)
      expect(json_body[:files][0][:encoding]).to eq("utf-8")
    end

    it 'handles binary file correctly' do
      workspace = Vulcan::Workspace.all[0]
      remote_manager.cp_file("/test-utils/example.bin","#{workspace.path}/output/example.bin")
      request = {
        file_names: ["example.bin"]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/read", request)
      expect(last_response.status).to eq(200)
      expect(json_body[:files][0][:encoding]).to eq("base64")
    end

    it 'handles json correctly' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      file_name = "test_file_name"
      json_content = '{"contents":[{"col":"gene1","def":["exactly",0,"below",0.02],"logic":null}]}'
      request = {
        files: [{
          filename: file_name,
          content: json_content
        }]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/write", request)
      expect(last_response.status).to eq(200)

      request = {
        file_names: ["test_file_name"]
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/read", request)
      expect(json_body[:files][0][:content]).to eq(json_content)
      expect(json_body[:files][0][:encoding]).to eq("utf-8")
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

  context 'get config' do
    before do
      setup_workspace
    end
    
    it 'gets the config of the workspace' do
      # Send a config for the first job
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
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      expect(last_response.status).to eq(200)

      # Get the config of the workspace
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config/#{config_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:id]).to_not be_nil
      expect(json_body[:workspace_id]).to_not be_nil
      expect(json_body[:path]).to_not be_nil
      expect(json_body[:hash]).to_not be_nil
      expect(json_body[:input_files]).to_not be_nil
      expect(json_body[:input_params]).to_not be_nil
      expect(json_body[:state]).to_not be_nil
    end

    it 'returns an error if the config does not exist' do
      auth_header(:editor)
      workspace_id = Vulcan::Workspace.all[0].id
      get("/api/v2/#{PROJECT}/workspace/#{workspace_id}/config/1234567890")
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("Config not found")
    end

  end

  context 'get state' do
    before do
      setup_workspace
    end

    it 'gets the state for a valid config' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      
      # Create a config
      request = {
        params: {
          count_bytes: false,
          count_chars: true
        },
        uiFilesSent: [],
        paramsChanged: []
      }
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      expect(last_response.status).to eq(200)
      
      # Store the state returned from save-config
      save_config_state = {
        files: json_body[:files],
        jobs: json_body[:jobs]
      }

      # Get the state for this config
      get("/api/v2/#{PROJECT}/config/#{config_id}/state")
      expect(last_response.status).to eq(200)
      
      # Verify the state contains expected files and jobs
      expect(json_body[:files][:planned]).to match_array(["output/count_poem.txt", "output/count_poem_2.txt"])
      expect(json_body[:jobs][:planned]).to match_array(["count"])
      
      # Verify that the state from get-state matches the state from save-config
      expect(json_body[:files]).to eq(save_config_state[:files])
      expect(json_body[:jobs]).to eq(save_config_state[:jobs])
    end

    it 'returns an error if the config does not exist' do
      auth_header(:editor)
      get("/api/v2/#{PROJECT}/config/1234567890/state")
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to include("does not exist")
    end

    it 'returns different state after we have run the config' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      
      # Create a config
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
      config_id = json_body[:config_id]
      
      # Store the initial state from save-config
      initial_state = {
        files: json_body[:files],
        jobs: json_body[:jobs]
      }
      
      # Run the workflow
      run_with_retry do
        post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      end
      expect(last_response.status).to eq(200)
      
      # Wait for the job to complete
      run_id = json_body[:run_id]
      check_jobs_status(["count"]) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
      end
      expect(last_response.status).to eq(200)

      # Get the state after running the workflow
      get("/api/v2/#{PROJECT}/config/#{config_id}/state")
      expect(last_response.status).to eq(200)
      
      # Verify that the state has changed after running the workflow
      # Files that were planned should now be completed
      current_state = json_body
      expect(current_state[:files][:completed]).to match_array(["output/count_poem.txt", "output/count_poem_2.txt"])
      expect(current_state[:jobs][:completed]).to match_array(["count"])
      
      # The completed files/jobs should be different from the initial planned state
      expect(current_state[:files][:completed]).to_not eq(initial_state[:files][:completed])
      expect(current_state[:jobs][:completed]).to_not eq(initial_state[:jobs][:completed])
      
      # The planned arrays should be empty since the job completed
      expect(current_state[:files][:planned]).to be_empty
      expect(current_state[:jobs][:planned]).to be_empty
    end
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

  context 'status checking' do

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
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/config", request)
      config_id = json_body[:config_id]
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{config_id}")
      expect(last_response.status).to eq(200)
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{json_body[:run_id]}")
      expect(last_response.status).to eq(200)
      expect(["NOT STARTED", "RUNNING"]).to include(json_body[:count])
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

  context 'download file' do

    before do
      setup_workspace
    end

    it 'returns the file for download' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_files_to_workspace(workspace.id)
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/download/poem.txt")
      expect(last_response.status).to eq(200)
      # Check headers
      expect(last_response.headers["content-type"]).to eq({file_type: "application/octet-stream", disposition: "attachment; filename=poem.txt"})
      expect(last_response.headers["content-length"]).to eq("142")
    end

    it 'returns a file nested in a directory' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      write_nested_file_to_workspace(workspace.id)
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/download/nested_dir/poem.txt")
      expect(last_response.status).to eq(200)
      # Check headers for nested file
      expect(last_response.headers["content-type"]).to eq({file_type: "application/octet-stream", disposition: "attachment; filename=nested_dir/poem.txt"})
      expect(last_response.headers["content-length"]).to eq("142")
    end

    it 'raises an error if the file does not exist' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/download/does_not_exist.txt") 
      expect(last_response.status).to eq(422)
    end
  end

  context 'stream file' do

    before do
      setup_workspace
    end

    it 'streams the file for download' do
      auth_header(:editor)

      workspace = Vulcan::Workspace.first
      write_files_to_workspace(workspace.id)

      get "/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/stream-download/poem.txt"

      # --- force the body to be consumed (two equivalent ways) ---
      streamed = +""
      last_response.body.each { |chunk| streamed << chunk }      # explicit
      # streamed = last_response.body                            # implicit join

      expect(last_response.status).to eq(200)

      # headers: rack sets canonical case
      expect(last_response.headers['Content-Type']).to eq('application/octet-stream')
      expect(last_response.headers['Content-Disposition']).to eq('attachment; filename=poem.txt')
      expect(last_response.headers['Cache-Control']).to eq('no-cache')
    end

    it 'raises an error if the file does not exist' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/file/stream-download/does_not_exist.txt") 
      expect(last_response.status).to eq(422)
    end
  end

  context 'cluster latency' do
    before do
      Vulcan.instance.remove_instance_variable('@latency') if Vulcan.instance.instance_variable_defined?('@latency')
    end

    it 'fails without auth' do
      get("/api/v2/cluster-latency")
      expect(last_response.status).to eq(401)
    end

    it 'returns SSH latency measurement' do
      auth_header(:editor)
      get("/api/v2/cluster-latency")
      expect(last_response.status).to eq(200)
      expect(json_body[:latency]).to be_a(Float)
    end

    it 'actually measures latency' do
      allow_any_instance_of(Vulcan::RemoteManager).to receive(:measure_latency).and_return(2)
      auth_header(:editor)
      get("/api/v2/cluster-latency")
      expect(last_response.status).to eq(200)
      expect(json_body[:latency]).to be_a(Float)
      expect(Vulcan.instance.instance_variable_get('@remote_manager')).to have_received(:measure_latency)
    end

    context 'with caching' do 
      it 'skips measuring latency if unnecessary' do
        manager = Vulcan::RemoteManager.new(Vulcan.instance.ssh_pool)
        allow(manager).to receive(:measure_latency).and_return(2)
        Vulcan.instance.instance_variable_set("@remote_manager", manager)
        Vulcan.instance.instance_variable_set("@latency", 5.times.map do { latency: 2, date: Time.now } end)

        auth_header(:editor)
        get("/api/v2/cluster-latency")
        expect(last_response.status).to eq(200)
        expect(manager).not_to have_received(:measure_latency)
      end

      it 'does not recomputes if forced and data is current' do
        manager = Vulcan::RemoteManager.new(Vulcan.instance.ssh_pool)
        allow(manager).to receive(:measure_latency).and_return(2)
        Vulcan.instance.instance_variable_set("@remote_manager", manager)
        Vulcan.instance.instance_variable_set("@latency", 5.times.map do { latency: 2, date: Time.now } end)

        auth_header(:editor)
        get("/api/v2/cluster-latency?recompute=true")
        expect(last_response.status).to eq(200)
        expect(manager).not_to have_received(:measure_latency)
      end

      it 'recomputes if forced and data is out of date' do
        manager = Vulcan::RemoteManager.new(Vulcan.instance.ssh_pool)
        allow(manager).to receive(:measure_latency).and_return(2)
        Vulcan.instance.instance_variable_set("@remote_manager", manager)
        Vulcan.instance.instance_variable_set("@latency", 5.times.map do { latency: 2, date: Time.now - Vulcan.instance.config(:latency_time) * 2 } end)

        auth_header(:editor)
        get("/api/v2/cluster-latency?recompute=true")
        expect(last_response.status).to eq(200)
        expect(manager).to have_received(:measure_latency)
      end
    end
  end

  context 'cancel workflow' do
    before do
      setup_workspace
    end

    def start_workflow(workspace)
      auth_header(:editor)
      write_files_to_workspace(workspace.id)
      
      # Start a workflow
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
      expect(last_response.status).to eq(200)
      run_id = json_body[:run_id]
      return run_id
    end

    xit 'successfully cancels a running workflow' do
      workspace = Vulcan::Workspace.all[0]
      run_id = start_workflow(workspace)

      # Verify the workflow is running
      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/running")
      expect(last_response.status).to eq(200)
      expect(json_body[:running]).to be_truthy

      get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")

      # Cancel the workflow
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/#{run_id}/cancel")
      expect(last_response.status).to eq(200)
      expect(json_body[:message]).to eq("Workflow canceled successfully")
      expect(json_body[:run_id]).to eq(run_id)

      # To guard against race conditions, we immediately capture the job statuses
      # is_running just checks that the lock file has been removed
      # which is a sign that snakemake has self cleaned up
      job_status_history = { count: [], arithmetic: [], checker: [] }

      retry_until(max_attempts: 5, base_delay: 5) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/run/#{run_id}")
        job_statuses = json_body
        
        # Append current status to history for each job
        job_statuses.each do |job_name, status|
          if job_status_history[job_name.to_sym]
            job_status_history[job_name.to_sym] << status
          end
        end
        
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/running")
        expect(last_response.status).to eq(200)
        json_body[:running] == false
      end
      expect(json_body[:running]).to be_falsey

      # Verify that at least one job was cancelled
      expect(last_response.status).to eq(200)
      cancelled_jobs = job_status_history.values.flatten.select { |status| status == "CANCELLED by 0" }
      expect(cancelled_jobs).not_to be_empty
    end

    it 'fails to cancel when no workflow is running' do
      workspace = Vulcan::Workspace.all[0]
      run_id = start_workflow(workspace)
      # Cancel the workflow
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/#{run_id}/cancel")
      expect(last_response.status).to eq(200)
      expect(json_body[:message]).to eq("Workflow canceled successfully")
      # Wait till its cancelled
      retry_until(max_attempts: 5, base_delay: 5) do
        get("/api/v2/#{PROJECT}/workspace/#{workspace.id}/running")
        expect(last_response.status).to eq(200)
        json_body[:running] == false
      end
      
      # Cancel the workflow again
      post("/api/v2/#{PROJECT}/workspace/#{workspace.id}/#{run_id}/cancel")
      expect(last_response.status).to eq(200)
      expect(json_body[:message]).to include("No workflow is currently running")
    end

  end

  context 'delete workspace' do
    before do
      setup_workspace
    end

    it 'successfully deletes a workspace when user is admin' do
      auth_header(:admin)
      workspace = Vulcan::Workspace.all[0]
      workspace_id = workspace.id
      workspace_path = workspace.path
      
      # Delete the workspace
      delete("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:message]).to eq("Workspace #{workspace.name} deleted successfully")
      expect(json_body[:workspace_id]).to eq(workspace_id.to_s)
      
      # Verify workspace directory is deleted
      expect(remote_manager.dir_exists?(workspace_path)).to be_falsey

      # Workspace, configs, and runs should be deleted
      expect(Vulcan::Workspace.first(id: workspace_id)).to be_nil
      expect(Vulcan::Config.where(workspace_id: workspace_id).count).to eq(0)
      expect(Vulcan::Run.where(workspace_id: workspace_id).count).to eq(0)
    end

    it 'fails to delete workspace when user is not admin' do
      auth_header(:editor)
      workspace = Vulcan::Workspace.all[0]
      workspace_id = workspace.id
      
      # Try to delete the workspace
      delete("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq("You are forbidden from performing this action.")
      
      # Verify workspace still exists
      expect(Vulcan::Workspace.first(id: workspace_id)).to_not be_nil
    end

    it 'fails to delete workspace when user is viewer' do
      auth_header(:viewer)
      workspace = Vulcan::Workspace.all[0]
      workspace_id = workspace.id
      
      # Try to delete the workspace
      delete("/api/v2/#{PROJECT}/workspace/#{workspace_id}")
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq("You are forbidden from performing this action.")
      
      # Verify workspace still exists
      expect(Vulcan::Workspace.first(id: workspace_id)).to_not be_nil
    end

    it 'fails to delete non-existent workspace' do
      auth_header(:admin)
      non_existent_id = 99999
      
      # Try to delete non-existent workspace
      delete("/api/v2/#{PROJECT}/workspace/#{non_existent_id}")
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to include("does not exist")
    end

  end
end
