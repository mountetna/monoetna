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
      get("/api/#{PROJECT}/workflows")

      expect(last_response.status).to eq(200)

      response = JSON.parse(last_response.body)
      expect(response['workflows'].first['inputs']).to eql({
          "someInt" => {
              "default" => 200,
              "format" => nil,
              "label" => "it is an int",
              "type" => "int",
          },
          "someIntWithoutDefault" => {
              "default" => nil,
              "format" => nil,
              "label" => nil,
              "type" => "int",
          },
      })

      expect(response['workflows'].first['outputs']).to eql({
          "the_result" => {
              "default" => nil,
              "format" => nil,
              "label" => nil,
              "outputSource" => "finalStep/sum",
              "type" => "int",
          }
      })

      expect(response['workflows'].first['steps']).to eql([
          [
              {
                  "in" => ["a", "b"],
                  "out" => ["sum"],
                  "name" => "firstAdd",
                  "run" => "scripts/add.cwl",
              },
              {
                  "in" => ["a", "b"],
                  "out" => ["sum"],
                  "name" => "finalStep",
                  "run" => "scripts/add.cwl",
              },
          ],
          [
              {
                  "in" => ["a", "b"],
                  "out" => ["sum"],
                  "name" => "firstAdd",
                  "run" => "scripts/add.cwl",
              },
              {
                  "in" => ["num"],
                  "out" => ["num"],
                  "name" => "pickANum",
                  "run" => "ui-queries/pick-a-number.cwl",
              },
              {
                  "in" => ["a", "b"],
                  "out" => ["sum"],
                  "name" => "finalStep",
                  "run" => "scripts/add.cwl",
              },
          ]
      ])
    end

    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/#{PROJECT}/workflows")

      expect(last_response.status).to eq(403)
    end
  end
end
