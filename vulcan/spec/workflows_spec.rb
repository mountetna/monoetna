require_relative '../lib/server/controllers/workflows_controller'

describe WorkflowsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  describe 'loading' do
    ::Dir.glob(::File.join(::File.dirname(__FILE__), '..', 'lib', 'server', 'workflows', '*.{yaml,yml,cwl}')).each do |yml|
      it "#{yml} works" do
        Etna::Cwl::Workflow.from_yaml_file(File.basename(yml), File.dirname(yml))
      end
    end
  end

  context '#fetch' do
    it 'does not list workflows the user cannot access' do
      auth_header(:viewer, additional: { perm: "v:not-a-thing"})
      get("/api/workflows")

      expect(last_response.status).to eq(200)
      response = JSON.parse(last_response.body)

      workflows = response['workflows'].map { |w| w['name'] }
      # This workflow didn't have any projects specified, so it still shows up
      # However, other workflows are not included.
      expect(workflows).to eql(["test_concurrent_workflow.cwl"])
    end

    it 'gets a list of workflows' do
      auth_header(:viewer)
      get("/api/workflows")
      expect(last_response.status).to eq(200)

      response = JSON.parse(last_response.body)
      workflow = response['workflows'].find { |w| w['name'] == 'test_workflow.cwl' }

      expect(workflow['inputs']).to eql({
          "someInt" => {
              "default" => 200,
              "format" => nil,
              "label" => "it is an int",
              "type" => "int",
              "doc" => "help tip"
          },
          "someIntWithoutDefault" => {
              "default" => nil,
              "format" => nil,
              "label" => nil,
              "type" => "int",
              "doc" => "another tip"
          },
      })

      # does not include the secret project
      expect(workflow['projects']).to eql(['labors'])

      expect(workflow['outputs']).to eql({
          "the_result" => {
              "default" => nil,
              "format" => nil,
              "label" => nil,
              "outputSource" => "finalStep/sum",
              "type" => "int",
          },
          "thumbnail" => {
              "default" => nil,
              "format" => "image/png",
              "label" => nil,
              "outputSource" => "finalStep/thumb.png",
              "type" => "File"
          }
      })

      expect(workflow['steps']).to eql([
          [
              {
                  "in" => [{"id"=>"a", "source"=>"someInt"},
                      {"id"=>"b", "source"=>"someIntWithoutDefault"}],
                  "doc" => nil,
                  "label"=>nil,
                  "out" => ["sum"],
                  "name" => "firstAdd",
                  "run" => "scripts/add.cwl",
              },
              {
                  "in" => [{"id"=>"num", "source"=>"firstAdd/sum"}],
                  "doc" => nil,
                  "label"=>nil,
                  "out" => ["num"],
                  "name" => "pickANum",
                  "run" => "ui-queries/pick-a-number.cwl",
              },
              {
                  "in" => [{"id"=>"a", "source"=>"firstAdd/sum"},
                      {"id"=>"b", "source"=>"pickANum/num"}],
                  "doc" => nil,
                  "label"=>nil,
                  "out" => ["sum", "thumb.png"],
                  "name" => "finalStep",
                  "run"=>"scripts/add.cwl"},
              {
                  "in" => [{"id"=>"a", "source"=>"finalStep/sum"},
                           {"id"=>"b", "source"=>"finalStep/thumb.png"}],
                  "doc" => nil,
                  "label"=>nil,
                  "name"=>"aPlot",
                  "out"=>[],
                  "run"=>"ui-outputs/plotly.cwl"}

          ],
      ])

      save_last_response_json('workflows-response', 'WorkflowsResponse')
    end

    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/workflows")

      expect(last_response.status).to eq(403)
    end

    it 'allows a user without vulcan flag' do
      auth_header(:no_flag)
      get("/api/workflows")

      expect(last_response.status).to eq(200)
    end
  end

  context '#fetch_for_project' do
    it 'does not list workflows the user cannot access' do
      auth_header(:viewer, additional: { perm: "v:not-a-thing"})
      get("/api/labors/workflows")

      expect(last_response.status).to eq(403)
    end

    it 'gets a list of workflows for the specific project' do
      auth_header(:viewer, additional: { perm: "v:labors"})
      get("/api/labors/workflows")
      expect(last_response.status).to eq(200)

      response = JSON.parse(last_response.body)
      workflows = response['workflows'].map { |w| w['name'] }

      expect(workflows).to match_array(['test_concurrent_workflow.cwl', 'test_workflow.cwl'])

      save_last_response_json('project-workflows-response', 'ProjectWorkflowsResponse')
    end

    it 'can filter project workflows by tag' do
      auth_header(:viewer, additional: { perm: "v:labors"})
      get("/api/labors/workflows?tag=demo")
      expect(last_response.status).to eq(200)

      response = JSON.parse(last_response.body)
      workflows = response['workflows'].map { |w| w['name'] }

      expect(workflows).to match_array(['test_workflow.cwl'])
    end

    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/labors/workflows")

      expect(last_response.status).to eq(403)
    end

    it 'allows a user without vulcan flag' do
      auth_header(:no_flag)
      get("/api/labors/workflows")

      expect(last_response.status).to eq(200)
    end
  end
end
