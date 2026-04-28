describe UpdateModelController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    stub_event_log
  end

  context "requests from non-superusers" do
    it "rejects the requests" do
      auth_header(:editor)

      json_post(:update_model, {
        project_name: "labors",
        actions: [
          action_name: "update_attribute",
          model_name: "monster",
          attribute_name: "name",
          description: "The monster's name",
        ]
      })

      expect(last_response.status).to eq(403)
    end
  end

  context "valid and authorized request" do
    before(:context) do
      ensure_labors_template_project
    end

    before do
      @original_attribute = Labors::Monster.attributes[:name].dup
    end

    after do
      attribute = Labors::Monster.attributes[:name]
      attribute.description = @original_attribute.description
      attribute.display_name = @original_attribute.display_name
    end

    it "returns the project template with all changes" do
      auth_header(:superuser)
      json_post(:update_model, {
        project_name: "labors",
        actions: [{
          action_name: "update_attribute",
          model_name: "monster",
          attribute_name: "name",
          description: "The monster's name",
          display_name: "NAME"
        }]
      })

      expect(last_response.status).to eq(200)

      response_json = JSON.parse(last_response.body)
      attribute_json = response_json["models"]["monster"]["template"]["attributes"]["name"]
      expect(attribute_json["description"]).to eq("The monster's name")
      expect(attribute_json["display_name"]).to eq("NAME")
    end

    it "returns model template links after a model metadata update" do
      auth_header(:superuser)
      json_post(:update_model, {
        project_name: "labors",
        actions: [{
          action_name: "set_model_template",
          model_name: "monster",
          template_project_name: "labors_template",
          template_model_name: "monster"
        }]
      })

      expect(last_response.status).to eq(200)

      response_json = JSON.parse(last_response.body)
      model_json = response_json["models"]["monster"]["template"]
      expect(model_json["template_project_name"]).to eq("labors_template")
      expect(model_json["template_model_name"]).to eq("monster")
    end
  end

  context "invalid action" do
    it "does not update attribute options with invalid attribute name" do
      auth_header(:superuser)
      json_post(:update_model, {
        project_name: "labors",
        actions: [{
          action_name: "update_attribute",
          model_name: "monster",
          attribute_name: "something_invalid",
          display_name: "NAME"
        }]
      })

      response_json = JSON.parse(last_response.body)
      expect(last_response.status).to eq(422)
      expect(response_json['errors'][0]['message']).to eq('Attribute does not exist')

      # No other changes are made when actions fail
      expect(Labors::Monster.attributes[:name].display_name).not_to eq("NAME")
    end
  end
end
