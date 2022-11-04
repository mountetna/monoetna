describe Magma::ReparentModelAction do
  let(:action) { Magma::ReparentModelAction.new("labors", action_params) }

  before(:each) { Timecop.freeze('2000-01-01') } # 946684800
  after(:each) { Timecop.return }

  describe "#perform" do
    context "for child/collection parent_link_type" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "sidekick",
          parent_model_name: "monster"
        }
      end

      before do
        @original_parent_attribute = Labors::Victim.attributes[:sidekick].dup
      end

      after do
        # Re-add parent relationship
        Labors::Victim.attributes[:sidekick] = @original_parent_attribute
      end

      it "removes the model link and attaches it to the new parent" do
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
          action_name: "reparent_model",
          model_name: "treatment",
          parent_model_name: "monster"
        }
      end

      before do
        @original_parent_attribute = Labors::Wound.attributes[:treatment].dup
      end

      after do
        # Re-add parent relationship
        Labors::Wound.attributes[:treatment] = @original_parent_attribute
      end

      it "removes the model link and attaches it to the new parent" do
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
          action_name: "reparent_model",
          model_name: "vegetation",
          parent_model_name: "victim"
        }
      end

      before do
        @original_parent_attribute = Labors::Habitat.attributes[:vegetation].dup
      end

      after do
        # Re-add parent relationship
        Labors::Habitat.attributes[:vegetation] = @original_parent_attribute
      end

      it "removes the model link and attaches it to the new parent" do
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
          action_name: "reparent_model",
          model_name: "rodent",
          parent_model_name: "monster"
        }
      end

      before do
        @original_parent_attribute = Labors::Habitat.attributes[:rodent].dup
      end

      after do
        # Re-add parent relationship
        Labors::Habitat.attributes[:rodent] = @original_parent_attribute
      end

      it "removes the model link and attaches it to the new parent" do
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
    let(:user) {
      Etna::User.new(
        AUTH_USERS[:superuser],
        'some.test.token'
      )
    }

    before(:each) do
      @project = create(:project, name: 'The Twelve Labors of Hercules')
    end

    context "when model has data" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "sidekick",
          parent_model_name: "monster",
          user: user
        }
      end

      before(:each) do
        lion = create(:labor, name: 'Nemean Lion', number: 1, project: @project)
        lion_monster = create(:monster, name: 'Nemean Lion', labor: lion)
        victim_1 = create(:victim, name: 'John Doe', monster: lion_monster)
        sidekick_1 = create(:sidekick, name: 'Jane Doe', victim: victim_1)
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
      end
    end

    context "when model has restricted data and user cannot see it" do
      let(:user) {
        Etna::User.new(
          AUTH_USERS[:editor],
          'some.test.token'
        )
      }

      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "sidekick",
          parent_model_name: "monster",
          user: user
        }
      end

      before(:each) do
        lion = create(:labor, name: 'Nemean Lion', number: 1, project: @project)
        lion_monster = create(:monster, name: 'Nemean Lion', labor: lion)
        victim_1 = create(:victim, name: 'John Doe', monster: lion_monster)
        sidekick_1 = create(:sidekick, name: 'Jane Doe', victim: victim_1, restricted: true)
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
      end
    end

    context "when child model has data" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "victim",
          parent_model_name: "habitat",
          user: user
        }
      end

      before(:each) do
        # a disconnect record, perhaps
        sidekick_1 = create(:sidekick, name: 'Jane Doe')
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
      end
    end

    context "when child model has restricted data and user cannot see it" do
      let(:user) {
        Etna::User.new(
          AUTH_USERS[:editor],
          'some.test.token'
        )
      }

      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "victim",
          parent_model_name: "habitat",
          user: user
        }
      end

      before(:each) do
        sidekick_1 = create(:sidekick, name: 'Jane Doe', restricted: true)
      end

      it "returns false and adds an error" do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
      end
    end

    context "when model does not exist" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "imaginary_model",
          parent_model_name: "monster"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Model does not exist.")
      end
    end

    context "when parent model does not exist" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "sidekick",
          parent_model_name: "imaginary_model"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Parent model does not exist.")
      end
    end

    context "when reparent causes a cyclical graph" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "victim",
          parent_model_name: "sidekick"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Action would lead to a cyclical graph, since parent_model_name is in the model_name tree.")
      end
    end

    context "when model same as parent_model" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "victim",
          parent_model_name: "victim"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Action would lead to a cyclical graph, since parent_model_name is in the model_name tree.")
      end
    end

    context "when parent_model doesn't change" do
      let(:action_params) do
        {
          action_name: "reparent_model",
          model_name: "victim",
          parent_model_name: "monster"
        }
      end

      it 'captures an error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("parent_model_name is already the parent of model_name.")
      end
    end
  end
end
