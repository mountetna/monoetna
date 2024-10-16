require_relative '../lib/magma'
require 'yaml'

describe Magma::Project do
  describe ".initialize" do
    it "loads attributes on the project's models from the database" do
      Magma.instance.db[:attributes].insert(
        project_name: "labors",
        model_name: "monster",
        attribute_name: "size",
        column_name: "size",
        type: "string",
        created_at: Time.now,
        updated_at: Time.now,
        display_name: "Monster Size",
        description: "How big is this monster?",
        unique: true
      )

      Magma.instance.db[:attributes].insert(
        project_name: "labors",
        model_name: "monster",
        attribute_name: "color",
        column_name: "color",
        type: "string",
        created_at: Time.now,
        updated_at: Time.now,
        description: "What color is it?",
      )

      project = Magma::Project.new(project_name: "labors")
      # Fetch and delete test attributes so they don't affect other tests
      size_attribute = Labors::Monster.attributes.delete(:size)
      color_attribute = Labors::Monster.attributes.delete(:color)

      expect(size_attribute.display_name).to eq("Monster Size")
      expect(size_attribute.description).to eq("How big is this monster?")
      expect(size_attribute.unique).to be(true)
      expect(color_attribute.description).to eq("What color is it?")
    end

    it "loads link models defined in the database" do
      Magma.instance.db[:attributes].insert(
        project_name: "labors",
        model_name: "monster",
        attribute_name: "alter_ego",
        column_name: "alter_ego",
        type: "link",
        created_at: Time.now,
        updated_at: Time.now,
        link_model_name: "monster"
      )

      project = Magma::Project.new(project_name: "labors")
      # Fetch and delete test attributes so they don't affect other tests
      attribute = Labors::Monster.attributes.delete(:alter_ego)

      expect(attribute.link_model).to eq(Labors::Monster)
    end

    it "raises an error when the database has attributes for a model that doesn't exist" do
      Magma.instance.db[:attributes].insert(
        project_name: "labors",
        model_name: "ghost",
        attribute_name: "name",
        column_name: "name",
        type: "string",
        created_at: Time.now,
        updated_at: Time.now,
        description: "There isn't a Labors::Ghost model"
      )

      expect { Magma::Project.new(project_name: "labors") }.to(
        raise_error(Magma::Project::AttributeLoadError)
      )
    end

    it "loads the project model's from the database" do
      Magma.instance.db[:models].insert(
        project_name: "movies",
        model_name: "hero",
        created_at: Time.now,
        updated_at: Time.now
      )

      Magma.instance.db[:models].insert(
          project_name: "movies",
          model_name: "status",
          created_at: Time.now,
          updated_at: Time.now
      )

      Magma.instance.db[:attributes].insert(
        project_name: "movies",
        model_name: "hero",
        attribute_name: "name",
        column_name: "name",
        type: "identifier",
        created_at: Time.now,
        updated_at: Time.now,
        description: "The hero's name"
      )

      project = Magma::Project.new(project_name: :movies)

      expect(project.models[:hero]).to eq(Movies::Hero)
      expect(Movies::Hero.attributes[:name].description).to eq("The hero's name")

      expect(project.models[:status]).to eq(Movies::Status)
    end

    it "loads project flags from the database" do
      create(:flag, project_name: "movies", flag_name: "is_pg", value: "true")
      create(:flag, project_name: "movies", flag_name: "is_in_theaters", value: "false")

      flags = Magma::Project.flags("movies")

      expect(flags['is_pg']).to eq('true')
      expect(flags['is_in_theaters']).to eq('false')

    end
  end

  describe '.ordered_models' do
    let(:project) { Magma.instance.magma_projects[:labors] }
    let(:model) { project.models[:monster] }

    it 'returns ordered models for labors project' do
      expect(project.ordered_models(model)).to match_array([Labors::Victim, Labors::Aspect, Labors::Sidekick, Labors::Wound, Labors::Treatment])
    end
  end
end
