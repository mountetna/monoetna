describe Session do
  describe "#as_json <-> from_json" do
    it "works" do
      session = Session.new_session_for("project_name", "workflow", "thekey", reference_figure_id: 1)
      expect(Session.from_json(JSON.parse(session.as_json.to_json)).as_json).to eql(session.as_json)
    end
  end

  describe "#workflow" do
    context "when reference id exists, loads workflow" do
      let(:figure) { create_figure_with_snapshot }

      it "from snapshot when snapshot exists" do
        session = Session.new_session_for("project_name", "test_workflow.cwl", "thekey", reference_figure_id: figure.id)
        expect(session.workflow.inputs.length).to eq(3)
        expect(session.workflow.inputs.map { |i| i.id }).to include("removedInt")
      end

      it "from yaml file when no snapshot" do
        figure.remove_existing_snapshot

        session = Session.new_session_for("project_name", "test_workflow.cwl", "thekey", reference_figure_id: figure.id)
        expect(session.workflow.inputs.length).to eq(2)
        expect(session.workflow.inputs.map { |i| i.id }).not_to include("removedInt")
      end
    end

    context "without reference id" do
      it "loads workflow from yaml file" do
        session = Session.new_session_for("project_name", "test_workflow.cwl", "thekey")
        expect(session.workflow.inputs.length).to eq(2)
        expect(session.workflow.inputs.map { |i| i.id }).not_to include("removedInt")
      end
    end
  end
end

describe SessionsController, e2e: true do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:mock_create_task_token) { true }
  let(:storage) { Vulcan::Storage.new }
  before(:each) do
    FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)
  end

  let(:session) do
    if body_json.include?(:key)
      Session.from_json(JSON.parse(body).merge({ "workflow_name" => workflow_name, "project_name" => project_name }))
    else
      Session.from_json(last_json_response["session"])
    end
  end

  let(:orchestration) { session.orchestration }
  let(:headers) do
    {
      "CONTENT_TYPE" => "application/json",
    }
  end
  let(:project_name) { PROJECT }
  let(:workflow_name) { "test_workflow.cwl" }
  let(:workflow_path) { "/api/#{project_name}/session/#{workflow_name}" }
  let(:body) { body_json.to_json }
  let(:body_json) do
    {
      key: key,
      inputs: inputs,
    }
  end
  let(:key) { "mykey" }
  let(:inputs) { {} }

  def make_request(postfix = "")
    post("#{workflow_path}#{postfix}", body, headers)
  end

  def last_json_response
    JSON.parse(last_response.body)
  end

  before(:each) { auth_header(:viewer) }

  before(:each) do
    if mock_create_task_token
      token = Base64.strict_encode64(AUTH_USERS[:viewer].to_json)
      allow_any_instance_of(SessionsController).to receive(:create_task_token).and_return(token)
    end
  end

  describe "without a saved workflow snapshot" do
    describe "creating a new session" do
      let(:body_json) { {} }
      describe "without viewer permission" do
        before(:each) { auth_header(:non_user) }

        it "returns forbidden" do
          make_request("/status")
          expect(last_response.status).to eql(403)
        end
      end

      it "creates a new empty session and returns it" do
        make_request("/status")
        expect(last_response.status).to eql(200)
        expect(last_json_response["session"]["project_name"]).to eql(project_name)
        expect(last_json_response["session"]["key"]).to_not be_empty
        expect(last_json_response["session"]["workflow_name"]).to eql(workflow_name)
        expect(last_json_response["session"]["inputs"]).to eql({})
        expect(last_json_response["status"].map { |a| a.map { |v| v["name"] } }).to match_array([
          ["firstAdd", "pickANum", "finalStep", "aPlot"],
        ])
        expect(last_json_response["status"]).to eql([
          [
            { "downloads" => nil, "name" => "firstAdd", "status" => "pending", "hash" => orchestration.build_target_for("firstAdd").cell_hash },
            { "downloads" => nil, "name" => "pickANum", "status" => "pending", "hash" => orchestration.build_target_for("pickANum").cell_hash },
            { "downloads" => nil, "name" => "finalStep", "status" => "pending", "hash" => orchestration.build_target_for("finalStep").cell_hash },
            { "downloads" => nil, "name" => "aPlot", "status" => "pending", "hash" => orchestration.build_target_for("aPlot").cell_hash },
          ],
        ])
        expect(last_json_response["outputs"]).to eql({ "downloads" => nil, "status" => "pending" })

        save_last_response_json("status-without-downloads", "SessionStatusResponse")
      end
    end

    describe "with guest permission" do
      before(:each) { auth_header(:guest) }

      it "creates a new empty session and returns it" do
        make_request("/status")
        expect(last_response.status).to eql(200)
        expect(last_json_response["session"]["project_name"]).to eql(project_name)
        expect(last_json_response["session"]["key"]).to_not be_empty
        expect(last_json_response["session"]["workflow_name"]).to eql(workflow_name)
        expect(last_json_response["session"]["inputs"]).to eql({})
        expect(last_json_response["status"].map { |a| a.map { |v| v["name"] } }).to match_array([
          ["firstAdd", "pickANum", "finalStep", "aPlot"],
        ])
        expect(last_json_response["status"]).to eql([
          [
            { "downloads" => nil, "name" => "firstAdd", "status" => "pending", "hash" => orchestration.build_target_for("firstAdd").cell_hash },
            { "downloads" => nil, "name" => "pickANum", "status" => "pending", "hash" => orchestration.build_target_for("pickANum").cell_hash },
            { "downloads" => nil, "name" => "finalStep", "status" => "pending", "hash" => orchestration.build_target_for("finalStep").cell_hash },
            { "downloads" => nil, "name" => "aPlot", "status" => "pending", "hash" => orchestration.build_target_for("aPlot").cell_hash },
          ],
        ])
        expect(last_json_response["outputs"]).to eql({ "downloads" => nil, "status" => "pending" })

        save_last_response_json("status-without-downloads", "SessionStatusResponse")
      end
    end
  end

  describe "adding new inputs" do
    before(:each) do
      inputs["someIntWithoutDefault"] = 123
      inputs["pickANum/num"] = 300
    end

    def check_url_for(url, storage_file)
      get(URI.parse(url).path)
      expect(last_response["X-Sendfile"]).to eql(storage_file.data_path(storage))
    end

    it "builds and makes available downloads to those outputs" do
      make_request
      expect(last_response.status).to eql(200)
      response = last_json_response

      expect(response["status"][0][0]["status"]).to eql("running")
      expect(response["session"]["inputs"]).to eql(inputs)
      orchestration.scheduler.join_all

      make_request("/status")
      expect(last_response.status).to eql(200)
      save_last_response_json("status-with-downloads", "SessionStatusResponse")
      response = last_json_response

      check_url_for(response["status"][0][0]["downloads"]["sum"], orchestration.build_target_for("firstAdd").build_outputs["sum"])
      check_url_for(response["status"][0][1]["downloads"]["num"], orchestration.build_target_for("pickANum").build_outputs["num"])
      check_url_for(response["status"][0][2]["downloads"]["sum"], orchestration.build_target_for("finalStep").build_outputs["sum"])

      check_url_for(response["outputs"]["downloads"]["the_result"],
                    orchestration.build_target_for(:primary_outputs).build_outputs["the_result"])

      expect(response["status"][0][3]["status"]).to eql("complete")
    end

    describe "task_token" do
      describe "e2e" do
        # Set this to a project you have in your dev enviroment (not test environment) when re-recording
        let(:project_name) { "ipi" }
        let(:mock_create_task_token) { false }

        it "does a thing" do
          allow(Vulcan.instance).to receive(:config).and_call_original
          allow(Vulcan.instance).to receive(:config).with(:janus).and_return({
            host: "https://janus.development.local",
          })

          expect_any_instance_of(Vulcan::AsynchronousScheduler).to receive(:schedule_more!).and_wrap_original do |m, opts|
            new_token = opts[:token]
            payload = JSON.parse(Base64.decode64(new_token.split(".")[1]))
            # Ensure that the token is always downgraded to viewer in this case.
            expect(payload["perm"]).to eql("v:#{project_name}")
            expect(payload["task"]).to eql(true)

            m.call(opts)
          end

          # When re-recording, provide a TOKEN into the environment to run against your local development.
          if (tok = ENV["TOKEN"])
            header("Authorization", "Etna #{tok}")
          else
            # Since the project name will differ from a real project, when running this without a local token,
            # we need an auth token that still has permission to the project
            auth_header(:viewer, additional: { perm: "v:#{project_name}" })
          end

          VCR.use_cassette("create_session.e2e") do
            make_request
            expect(last_response.status).to eql(200)
          end
        end
      end
    end
  end

  describe "with a saved workflow snapshot" do
    describe "adding inputs from previous workflow version", use_transactional_fixtures: false do
      let(:figure) { create_figure_with_snapshot }

      let(:body_json) do
        {
          key: key,
          inputs: inputs,
          reference_figure_id: figure.id,
        }
      end

      before(:each) do
        inputs["someIntWithoutDefault"] = 123
        inputs["pickANum/num"] = 300
        inputs["removedInt"] = 42
      end

      def check_url_for(url, storage_file)
        get(URI.parse(url).path)
        expect(last_response["X-Sendfile"]).to eql(storage_file.data_path(storage))
      end

      it "builds and makes available downloads to those outputs" do
        make_request
        expect(last_response.status).to eql(200)
        response = last_json_response

        expect(response["status"][0][0]["status"]).to eql("running")
        expect(response["session"]["inputs"]).to eql(inputs)
        orchestration.scheduler.join_all

        make_request("/status")
        expect(last_response.status).to eql(200)
        save_last_response_json("status-with-downloads", "SessionStatusResponse")
        response = last_json_response

        check_url_for(response["status"][0][0]["downloads"]["sum"], orchestration.build_target_for("firstAdd").build_outputs["sum"])
        check_url_for(response["status"][0][1]["downloads"]["sum"], orchestration.build_target_for("secondAdd").build_outputs["sum"])
        check_url_for(response["status"][0][2]["downloads"]["num"], orchestration.build_target_for("pickANum").build_outputs["num"])
        check_url_for(response["status"][0][3]["downloads"]["sum"], orchestration.build_target_for("finalStep").build_outputs["sum"])

        check_url_for(response["outputs"]["downloads"]["the_result"],
                      orchestration.build_target_for(:primary_outputs).build_outputs["the_result"])

        expect(response["status"][0][4]["status"]).to eql("complete")
      end
    end
  end

  describe "handling errors" do
    before(:each) do
      inputs["someIntWithoutDefault"] = "abc"
      inputs["pickANum/num"] = "xyz"
    end

    def check_url_for(url, storage_file)
      get(URI.parse(url).path)
      expect(last_response["X-Sendfile"]).to eql(storage_file.data_path(storage))
    end

    it "reports error status and message" do
      make_request
      expect(last_response.status).to eql(200)
      orchestration.scheduler.join_all

      make_request("/status")
      expect(last_response.status).to eql(200)
      response = last_json_response

      expect(response["session"]["inputs"]).to eql(inputs)
      expect(response["status"].first[0]["status"]).to eq("error")
      expect(response["status"].first[2]["status"]).to eq("pending") # Can't run finalStep since firstStep has an error

      save_last_response_json("status-with-error", "SessionStatusResponse")
    end
  end
end
