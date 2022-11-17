describe Magma::RemoveModelAction do
  let(:action) { Magma::RemoveModelAction.new("labors", action_params) }

  before(:each) { Timecop.freeze('2000-01-01') } # 946684800
  after(:each) { Timecop.return }

  describe "#perform" do
    context "for child/collection parent_link_type" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "sidekick"
        }
      end

      before do
        @original_parent_attribute = Labors::Victim.attributes[:sidekick].dup
      end

      after do
        # Re-add parent relationship
        Labors::Victim.attributes[:sidekick] = @original_parent_attribute
      end

      it "removes the model link and detaches it from the graph" do
        project = Magma.instance.get_project(:labors)

        expect(project.detached_models.length).to eq(0)

        expect(action.perform).to eq(true)

        expect(project.detached_models.length).to eq(1)
        expect(project.detached_models).to match_array([Labors::Sidekick])
        expect(Labors::Victim.attributes[:sidekick]).to be_nil
      end
    end

    context "for table parent_link_type" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "treatment"
        }
      end

      before do
        @original_parent_attribute = Labors::Wound.attributes[:treatment].dup
      end

      after do
        # Re-add parent relationship
        Labors::Wound.attributes[:treatment] = @original_parent_attribute
      end

      it "removes the model link and detaches it from the graph" do
        project = Magma.instance.get_project(:labors)

        expect(project.detached_models.length).to eq(0)

        expect(action.perform).to eq(true)

        expect(project.detached_models.length).to eq(1)
        expect(project.detached_models).to match_array([Labors::Treatment])
        expect(Labors::Wound.attributes[:treatment]).to be_nil
      end
    end

    context "when has additional link up" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "vegetation"
        }
      end

      before do
        @original_parent_attribute = Labors::Habitat.attributes[:vegetation].dup
      end

      after do
        # Re-add parent relationship
        Labors::Habitat.attributes[:vegetation] = @original_parent_attribute
      end

      it "removes the model link and detaches it from the graph" do
        project = Magma.instance.get_project(:labors)

        expect(project.detached_models.length).to eq(0)

        expect(action.perform).to eq(true)

        expect(project.detached_models.length).to eq(1)
        expect(project.detached_models).to match_array([Labors::Vegetation])
        expect(Labors::Habitat.attributes[:vegetation]).to be_nil
      end
    end

    context "when has self-referential link" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "rodent"
        }
      end

      before do
        @original_parent_attribute = Labors::Habitat.attributes[:rodent].dup
      end

      after do
        # Re-add parent relationship
        Labors::Habitat.attributes[:rodent] = @original_parent_attribute
      end

      it "removes the model link and detaches it from the graph" do
        project = Magma.instance.get_project(:labors)

        expect(project.detached_models.length).to eq(0)

        expect(action.perform).to eq(true)

        expect(project.detached_models.length).to eq(1)
        expect(project.detached_models).to match_array([Labors::Rodent])
        expect(Labors::Habitat.attributes[:rodent]).to be_nil
      end
    end
  end

  describe "#validate" do
    context "when model has child" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "labor",
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model is not a leaf model.")
      end
    end

    context "when model has collection" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "habitat"
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model is not a leaf model.")
      end
    end

    context "when model has table" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "wound"
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model is not a leaf model.")
      end
    end

    context "when has additional link down" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "monster"
        }
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.map { |e| e[:message] }.include?("Model has a link to another model. Use \"remove_link\" on that attribute, first.")).to eq(true)
      end
    end

    context "when model does not exist" do
      let(:action_params) do
        {
          action_name: "remove_model",
          model_name: "imaginary_model"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model does not exist.")
      end
    end
  end
end
