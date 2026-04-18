describe Magma::SetModelTemplateAction do
  let(:action) { Magma::SetModelTemplateAction.new("labors", action_params) }

  describe "#perform" do
    let(:action_params) do
      {
        action_name: "set_model_template",
        model_name: "monster",
        template_project_name: "labors_template",
        template_model_name: "monster",
      }
    end

    it "sets the template model reference" do
      expect(action.perform).to eq(true)

      model_row = Magma.instance.db[:models].where(
        project_name: "labors",
        model_name: "monster"
      ).first

      expect(model_row[:template_project_name]).to eq("labors_template")
      expect(model_row[:template_model_name]).to eq("monster")
      expect(Labors::Monster.template_project_name).to eq("labors_template")
      expect(Labors::Monster.template_model_name).to eq(:monster)
    end

    it "clears the template model reference" do
      expect(action.perform).to eq(true)

      clear_action = Magma::SetModelTemplateAction.new("labors", {
        action_name: "set_model_template",
        model_name: "monster",
        clear_template: true,
      })

      expect(clear_action.perform).to eq(true)

      model_row = Magma.instance.db[:models].where(
        project_name: "labors",
        model_name: "monster"
      ).first

      expect(model_row[:template_project_name]).to eq(nil)
      expect(model_row[:template_model_name]).to eq(nil)
      expect(Labors::Monster.template_project_name).to eq(nil)
      expect(Labors::Monster.template_model_name).to eq(nil)
    end
  end

  describe "#validate" do
    context "when the model does not exist" do
      let(:action_params) do
        {
          action_name: "set_model_template",
          model_name: "iliad",
          template_project_name: "labors_template",
          template_model_name: "monster",
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model does not exist.")
      end
    end

    context "when template_model_name is missing" do
      let(:action_params) do
        {
          action_name: "set_model_template",
          model_name: "monster",
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Must include :template_model_name parameter")
      end
    end

    context "when clearing the template link" do
      let(:action_params) do
        {
          action_name: "set_model_template",
          model_name: "monster",
          clear_template: true,
        }
      end

      it "allows the action" do
        expect(action.validate).to eq(true)
      end
    end

    context "when the template model does not exist" do
      let(:action_params) do
        {
          action_name: "set_model_template",
          model_name: "monster",
          template_project_name: "labors_template",
          template_model_name: "ghost",
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Template model does not exist.")
      end
    end

    context "when the model points to itself" do
      let(:action_params) do
        {
          action_name: "set_model_template",
          model_name: "monster",
          template_project_name: "labors",
          template_model_name: "monster",
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model cannot point to itself as its template")
      end
    end
  end
end
