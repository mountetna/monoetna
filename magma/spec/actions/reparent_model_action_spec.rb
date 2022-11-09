# describe Magma::ReparentModelAction do
#   let(:action) { Magma::ReparentModelAction.new("labors", action_params) }
#   let(:revert_action) { Magma::ReparentModelAction.new("labors", revert_params) }

#   before(:each) { Timecop.freeze('2000-01-01') } # 946684800
#   after(:each) { Timecop.return }

#   describe "#perform" do
#     let(:project) { Magma.instance.get_project("labors") }
#     context "for child/collection parent_link_type" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "sidekick",
#           parent_model_name: "monster"
#         }
#       end
#       let(:revert_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "sidekick",
#           parent_model_name: "victim"
#         }
#       end

#       after do
#         revert_action.perform
#       end

#       it "removes the model link and attaches it to the new parent" do
#         expect(Labors::Monster.attributes[:sidekick]).to be_nil
#         expect(Labors::Sidekick.attributes[:victim]).not_to be_nil
#         expect(Labors::Sidekick.attributes[:monster]).to be_nil
#         expect(Labors::Victim.attributes[:sidekick]).not_to be_nil

#         expect(action.perform).to eq(true)

#         expect(Labors::Monster.attributes[:sidekick]).not_to be_nil
#         expect(Labors::Monster.attributes[:sidekick].class).to eq(Magma::CollectionAttribute)
#         expect(Labors::Sidekick.attributes[:victim]).to be_nil
#         expect(Labors::Sidekick.attributes[:monster]).not_to be_nil
#         expect(Labors::Sidekick.attributes[:monster].class).to eq(Magma::ParentAttribute)
#         expect(Labors::Victim.attributes[:sidekick]).to be_nil
#       end
#     end

#     context "for table parent_link_type" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "treatment",
#           parent_model_name: "monster"
#         }
#       end

#       let(:revert_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "treatment",
#           parent_model_name: "wound"
#         }
#       end

#       after do
#         revert_action.perform
#       end

#       it "removes the model link and attaches it to the new parent" do
#         expect(Labors::Monster.attributes[:treatment]).to be_nil
#         expect(Labors::Wound.attributes[:treatment]).not_to be_nil
#         expect(Labors::Treatment.attributes[:monster]).to be_nil
#         expect(Labors::Treatment.attributes[:wound]).not_to be_nil

#         expect(action.perform).to eq(true)

#         expect(Labors::Monster.attributes[:treatment]).not_to be_nil
#         expect(Labors::Monster.attributes[:treatment].class).to eq(Magma::TableAttribute)
#         expect(Labors::Wound.attributes[:treatment]).to be_nil
#         expect(Labors::Treatment.attributes[:monster]).not_to be_nil
#         expect(Labors::Treatment.attributes[:monster].class).to eq(Magma::ParentAttribute)
#         expect(Labors::Treatment.attributes[:wound]).to be_nil
#       end
#     end

#     context "when has additional link up" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "vegetation",
#           parent_model_name: "victim"
#         }
#       end

#       let(:revert_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "vegetation",
#           parent_model_name: "habitat"
#         }
#       end

#       after do
#         revert_action.perform
#       end

#       it "removes the model link and attaches it to the new parent" do
#         expect(Labors::Victim.attributes[:vegetation]).to be_nil
#         expect(Labors::Habitat.attributes[:vegetation]).not_to be_nil
#         expect(Labors::Vegetation.attributes[:victim]).to be_nil
#         expect(Labors::Vegetation.attributes[:habitat]).not_to be_nil

#         expect(action.perform).to eq(true)

#         expect(Labors::Victim.attributes[:vegetation]).not_to be_nil
#         expect(Labors::Habitat.attributes[:vegetation]).to be_nil
#         expect(Labors::Vegetation.attributes[:victim]).not_to be_nil
#         expect(Labors::Vegetation.attributes[:habitat]).to be_nil
#       end
#     end

#     context "when has self-referential link" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "rodent",
#           parent_model_name: "victim"
#         }
#       end

#       let(:revert_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "rodent",
#           parent_model_name: "habitat"
#         }
#       end

#       after do
#         revert_action.perform
#       end

#       it "removes the model link and attaches it to the new parent" do
#         expect(Labors::Victim.attributes[:rodent]).to be_nil
#         expect(Labors::Habitat.attributes[:rodent]).not_to be_nil
#         expect(Labors::Rodent.attributes[:habitat]).not_to be_nil

#         expect(action.perform).to eq(true)

#         expect(Labors::Victim.attributes[:rodent]).not_to be_nil
#         expect(Labors::Habitat.attributes[:rodent]).to be_nil
#         expect(Labors::Rodent.attributes[:habitat]).to be_nil
#       end
#     end
#   end

#   describe "#validate" do
#     let(:user) {
#       Etna::User.new(
#         AUTH_USERS[:superuser],
#         'some.test.token'
#       )
#     }

#     before(:each) do
#       @project = create(:project, name: 'The Twelve Labors of Hercules')
#     end

#     context "when model has data" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "sidekick",
#           parent_model_name: "monster",
#           user: user
#         }
#       end

#       before(:each) do
#         lion = create(:labor, name: 'Nemean Lion', number: 1, project: @project)
#         lion_monster = create(:monster, name: 'Nemean Lion', labor: lion)
#         victim_1 = create(:victim, name: 'John Doe', monster: lion_monster)
#         sidekick_1 = create(:sidekick, name: 'Jane Doe', victim: victim_1)
#       end

#       it "returns false and adds an error" do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
#       end
#     end

#     context "when model has restricted data and user cannot see it" do
#       let(:user) {
#         Etna::User.new(
#           AUTH_USERS[:editor],
#           'some.test.token'
#         )
#       }

#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "sidekick",
#           parent_model_name: "monster",
#           user: user
#         }
#       end

#       before(:each) do
#         lion = create(:labor, name: 'Nemean Lion', number: 1, project: @project)
#         lion_monster = create(:monster, name: 'Nemean Lion', labor: lion)
#         victim_1 = create(:victim, name: 'John Doe', monster: lion_monster)
#         sidekick_1 = create(:sidekick, name: 'Jane Doe', victim: victim_1, restricted: true)
#       end

#       it "returns false and adds an error" do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
#       end
#     end

#     context "when child model has data" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "victim",
#           parent_model_name: "habitat",
#           user: user
#         }
#       end

#       before(:each) do
#         # a disconnect record, perhaps
#         sidekick_1 = create(:sidekick, name: 'Jane Doe')
#       end

#       it "returns false and adds an error" do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
#       end
#     end

#     context "when child model has restricted data and user cannot see it" do
#       let(:user) {
#         Etna::User.new(
#           AUTH_USERS[:editor],
#           'some.test.token'
#         )
#       }

#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "victim",
#           parent_model_name: "habitat",
#           user: user
#         }
#       end

#       before(:each) do
#         sidekick_1 = create(:sidekick, name: 'Jane Doe', restricted: true)
#       end

#       it "returns false and adds an error" do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Cannot reparent a model with data records in its tree (including in child models).")
#       end
#     end

#     context "when model does not exist" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "imaginary_model",
#           parent_model_name: "monster"
#         }
#       end

#       it 'captures an error' do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Model does not exist.")
#       end
#     end

#     context "when parent model does not exist" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "sidekick",
#           parent_model_name: "imaginary_model"
#         }
#       end

#       it 'captures an error' do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Parent model does not exist.")
#       end
#     end

#     context "when reparent causes a cyclical graph" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "victim",
#           parent_model_name: "sidekick"
#         }
#       end

#       it 'captures an error' do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Action would lead to a cyclical graph, since parent_model_name is in the model_name tree.")
#       end
#     end

#     context "when model same as parent_model" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "victim",
#           parent_model_name: "victim"
#         }
#       end

#       it 'captures an error' do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("Action would lead to a cyclical graph, since parent_model_name is in the model_name tree.")
#       end
#     end

#     context "when parent_model doesn't change" do
#       let(:action_params) do
#         {
#           action_name: "reparent_model",
#           model_name: "victim",
#           parent_model_name: "monster"
#         }
#       end

#       it 'captures an error' do
#         expect(action.validate).to eq(false)
#         expect(action.errors.first[:message]).to eq("parent_model_name is already the parent of model_name.")
#       end
#     end
#   end
# end
