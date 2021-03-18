describe Session do
  describe '#as_json <-> from_json' do
    it 'works' do
      session = Session.new_session_for('project_name', 'workflow', 'thekey', {[:primary_inputs, "b"] => {json_payload: "a"}, ["primary_inputs", "c"] => {json_payload: "d"}})
      expect(Session.from_json(JSON.parse(session.as_json.to_json)).as_json).to eql(session.as_json)
    end
  end
end

describe SessionsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  let(:storage) { Vulcan::Storage.new }
  before(:each) do
    FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)
  end

  let(:session) { Session.from_json(JSON.parse(body).merge({ "workflow_name" => workflow_name, "project_name" => project_name })) }
  let(:orchestration) { session.orchestration }
  let(:headers) do
    {
        'CONTENT_TYPE' => 'application/json'
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
  let(:key) { 'mykey' }
  let(:inputs) { {} }

  def make_request(postfix = "")
    post("#{workflow_path}#{postfix}", body, headers)
  end

  def last_json_response
    JSON.parse(last_response.body)
  end

  before(:each) { auth_header(:viewer) }
  describe 'creating a new session' do
    let(:body_json) { {} }
    describe 'without viewer permission' do
      before(:each) do
        auth_header(:non_user)
      end

      it "returns forbidden" do
        make_request
        expect(last_response.status).to eql(403)
      end
    end

    it 'creates a new empty session and returns it' do
      make_request("/status")
      expect(last_response.status).to eql(200)
      expect(last_json_response['session']['project_name']).to eql(project_name)
      expect(last_json_response['session']['key']).to_not be_empty
      expect(last_json_response['session']['workflow_name']).to eql(workflow_name)
      expect(last_json_response['session']['inputs']).to eql({})
      expect(last_json_response['status'].map { |a| a.map { |v| v['name'] }}).to match_array([
          ['firstAdd', 'pickANum', 'finalStep', 'aPlot'],
      ])
      expect(last_json_response['status']).to match_array([
          [
              {'downloads' => nil, 'name' => 'firstAdd',  'status' => 'pending'},
              {'downloads' => nil, 'name' => 'pickANum',  'status' => 'pending'},
              {'downloads' => nil, 'name' => 'finalStep', 'status' => 'pending'},
              {'downloads' => nil, 'name' => 'aPlot',     'status' => 'pending'},
          ],
      ])
      expect(last_json_response['outputs']).to eql({'downloads' => nil, 'status' => 'pending'})
    end
  end

  describe 'adding new inputs' do
    before(:each) do
      inputs["someIntWithoutDefault"] = 123
      inputs["pickANum/num"] = 300
    end

    def check_url_for(url, storage_file)
      get(URI.parse(url).path)
      expect(last_response['X-Sendfile']).to eql(storage_file.data_path(storage))
    end

    it 'builds and makes available downloads to those outputs' do
      make_request
      expect(last_response.status).to eql(200)
      response = last_json_response

      expect(response['session']['inputs']).to eql(inputs)
      orchestration.scheduler.join_all

      make_request("/status")
      expect(last_response.status).to eql(200)
      response = last_json_response

      check_url_for(response['status'][0][0]['downloads']['sum'], orchestration.build_target_for('firstAdd').build_outputs['sum'])
      check_url_for(response['status'][0][1]['downloads']['num'], orchestration.build_target_for('pickANum').build_outputs['num'])
      check_url_for(response['status'][0][2]['downloads']['sum'], orchestration.build_target_for('finalStep').build_outputs['sum'])

      check_url_for(response['outputs']['downloads']['the_result'],
          orchestration.build_target_for(:primary_outputs).build_outputs['the_result'])

      expect(response['status'][0][3]['status']).to eql('complete')
    end
  end

  describe 'handling errors' do
    before(:each) do
      inputs["someIntWithoutDefault"] = 'abc'
      inputs["pickANum/num"] = 'xyz'
    end

    def check_url_for(url, storage_file)
      get(URI.parse(url).path)
      expect(last_response['X-Sendfile']).to eql(storage_file.data_path(storage))
    end

    it 'reports error status and message' do
      make_request
      expect(last_response.status).to eql(200)
      orchestration.scheduler.join_all

      make_request("/status")
      expect(last_response.status).to eql(200)
      response = last_json_response

      expect(response['session']['inputs']).to eql(inputs)
      expect(response['status'].first[0]['status']).to eq('error')
      expect(response['status'].first[2]['status']).to eq('pending') # Can't run finalStep since firstStep has an error
    end
  end
end
