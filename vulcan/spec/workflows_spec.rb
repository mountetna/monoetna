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

      expect(workflow['outputs']).to eql({
          "the_result" => {
              "default" => nil,
              "format" => nil,
              "label" => nil,
              "outputSource" => "finalStep/sum",
              "type" => "int",
          }
      })

      expect(workflow['dependencies_of_outputs']).to eql({
          "finalStep/sum" => [],
          "firstAdd/sum" => ["finalStep/sum", "pickANum/num"],
          "pickANum/num" => ["finalStep/sum"],
          "someInt" => ["firstAdd/sum", "finalStep/sum", "pickANum/num"],
          "someIntWithoutDefault" => ["firstAdd/sum", "finalStep/sum", "pickANum/num"],
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
                  "out" => ["sum"],
                  "name" => "finalStep",
                  "run"=>"scripts/add.cwl"},
              {
                  "in" => [{"id"=>"a", "source"=>"finalStep/sum"}],
                  "doc" => nil,
                  "label"=>nil,
                  "name"=>"aPlot",
                  "out"=>[],
                  "run"=>"ui-outputs/plotter.cwl"}
          ],
      ])

      save_last_response_json('workflows-response', 'WorkflowsResponse')
    end

    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/workflows")

      expect(last_response.status).to eq(403)
    end
  end
end
