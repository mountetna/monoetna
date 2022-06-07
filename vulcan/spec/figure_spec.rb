require_relative "../lib/server/controllers/figure_controller"

describe Vulcan::Figure do
  before(:each) do
    clear_store
  end

  it "returns a list of thumbnails" do
    store(
      "abee47d3ee8ba11e3fc2706d8d258e9e586465b4",
      "thumb.png",
      "thumbnail"
    )
    figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")

    expect(figure.thumbnails(storage: Vulcan::Storage.new)).to eq(["https://vulcan.test/api/labors/data/abee47d3ee8ba11e3fc2706d8d258e9e586465b4/thumb.png"])
  end

  it "ignores unbuilt thumbnails" do
    figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")

    expect(figure.thumbnails(storage: Vulcan::Storage.new)).to eq([])
  end
end

describe FigureController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#fetch" do
    it "returns a list of figures" do
      below_editor_roles.each do |role|
        auth_header(role)
        figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow")
        get("/api/labors/figures")

        expect(last_response.status).to eq(200)
        expect(json_body[:figures].first[:title]).to eql(figure.title)
      end
    end
  end

  context "#get" do
    it "returns a figure by id" do
      below_editor_roles.each do |role|
        auth_header(role)
        figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")
        get("/api/labors/figure/#{figure.figure_id}")

        expect(last_response.status).to eq(200)
        expect(json_body[:title]).to eql(figure.title)
        expect(json_body[:workflow_snapshot]).not_to eq(nil)
        expect(json_body[:workflow_snapshot][:steps][0][0].key?(:name)).to eq(true)
      end
    end


    it "snapshot includes metadata" do
      below_editor_roles.each do |role|
        auth_header(role)
        figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")
        get("/api/labors/figure/#{figure.figure_id}")

        expect(last_response.status).to eq(200)
        expect(json_body[:workflow_snapshot][:projects]).to eq(["labors", "secret_project"])
      end
    end
  end

  context "#revisions" do
    it "returns a figure's revisions'" do
      auth_header(:viewer)
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl", figure_id: 1)
      figure_v2 = create_figure(
        title: "Lion of Nemea, redux",
        workflow_name: "test_workflow.cwl",
        figure_id: figure.figure_id,
      )
      get("/api/labors/figure/#{figure.figure_id}/revisions")

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(2)
      expect(json_body.map { |r| r[:title] }).to match_array(["Lion of Nemea", "Lion of Nemea, redux"])
      expect(json_body.first.keys).to include(:id, :workflow_snapshot)
    end
  end

  context "#create" do
    it "creates a new figure" do
      expect(Vulcan::WorkflowSnapshot.count).to eq(0)
      below_editor_roles.each do |role|
        auth_header(role)
        contents = {
          title: "Lion of Nemea",
          workflow_name: "test_workflow.cwl",
          inputs: { a: "b" },
        }
        post("/api/labors/figure/create", contents)

        expect(last_response.status).to eq(200)
        expect(json_body).to include(
          contents.merge(
            project_name: "labors",
          )
        )
      end
    end

    it "creates a new public figure with tags" do
      auth_header(:viewer)
      contents = {
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        inputs: { a: "b" },
        tags: ["public"],
      }

      expect(Vulcan::Figure.count).to eq(0)
      expect(Vulcan::WorkflowSnapshot.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors",
        )
      )

      expect(json_body.keys).to include(:dependencies)
      expect(Vulcan::WorkflowSnapshot.count).to eq(1)
      expect(Vulcan::Figure.count).to eq(1)
    end

    it "refuses to create a public figure with tags for guest" do
      auth_header(:guest)
      contents = {
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        inputs: { a: "b" },
        tags: ["public"],
      }

      expect(Vulcan::Figure.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(422)

      expect(Vulcan::Figure.count).to eq(0)
    end

    it "guest can create a figure with non-public tags" do
      expect(Vulcan::Figure.count).to eq(0)
      auth_header(:guest)
      contents = {
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        inputs: { a: "b" },
        tags: ["mine"],
      }

      expect(Vulcan::Figure.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors",
        )
      )

      expect(json_body.keys).to include(:dependencies, :workflow_snapshot)
      expect(json_body[:workflow_snapshot]).not_to eq(nil)

      expect(Vulcan::Figure.count).to eq(1)
    end
  end

  context "#update" do
    it "updates an existing figure" do
      figure = create_figure(
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        dependencies: {
          something: "sha:abc",
        },
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)
      original_title = "Hercules Fighting the Nemean Lion"
      contents = { title: original_title, comment: "Better title" }

      below_editor_roles.each_with_index do |role, index|
        auth_header(role)
        contents[:title] = "#{original_title} #{index}"
        post("/api/labors/figure/#{figure.figure_id}/update", contents)

        expect(last_response.status).to eq(200)
        expect(json_body).to include(
          figure_id: 1,
          inputs: {},
          project_name: "labors",
          title: "Hercules Fighting the Nemean Lion #{index}",
          workflow_name: "test_workflow.cwl",

          # Updating maintains the previous figure dependencies
          dependencies: {
            something: "sha:abc",
          },
        )

        # there are more figures
        expect(Vulcan::Figure.count).to eq(2 + index)
        expect(Vulcan::WorkflowSnapshot.count).to eq(2 + index)
        expect(Vulcan::WorkflowSnapshot.first.as_steps_json_w_metadata).to eq(
          Vulcan::WorkflowSnapshot.last.as_steps_json_w_metadata
        )
      end
    end

    it "throws exception for unknown figure" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")
      auth_header(:viewer)
      contents = { title: "Hercules Fighting the Nemean Lion", comment: "Better title" }
      post("/api/labors/figure/999999999/update", contents)

      expect(last_response.status).to eq(404)
    end

    it "updates an existing figure with tags" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")
      auth_header(:viewer)
      contents = { tags: ["private"], comment: "Make private" }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        tags: ["private"],
        workflow_name: "test_workflow.cwl",
      )
    end

    it "can update dependencies" do
      figure = create_figure(
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        dependencies: {
          something: "sha:abc",
        },
      )

      auth_header(:viewer)
      contents = {
        title: "Hercules Fighting the Nemean Lion",
        comment: "Better title",
        update_dependencies: true,
      }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        title: "Hercules Fighting the Nemean Lion",
        workflow_name: "test_workflow.cwl",
      )

      # Using the update_dependencies flag should update dependencies
      expect(json_body[:dependencies]).not_to eq({
        something: "sha:abc",
      })

      # CI and local don't necessarily have production-tagged images, so
      #   the SHA's will just come back as empty strings.
      # But this is what should happen.
      # expect(json_body[:dependencies].values.all? { |d| d =~ /^sha256:.*/ }).to eq(true)
    end

    def setup_workflow_mocks
      allow(Etna::Cwl::Workflow).to receive(:metadata).and_return({
        "authors": ["a new author"],
        "projects": ["another project"],
      })
      allow(Etna::Cwl::Workflow).to receive(:raw_yaml_from_file).and_return(
        YAML.safe_load("---\ncwlVersion: v1.1\nclass: Workflow\ninputs: []\noutputs: []\nsteps:\n  dummyStep:\n    run: scripts/something.cwl\n    in: []\n    out: []\nname: test_workflow.cwl")
      )
      allow_any_instance_of(Vulcan::DependencyManager).to receive(:dependency_shas).and_return({
        "vulcan": "sha256:123",
      })
    end

    it "captures any workflow updates during figure update, if dependencies updated" do
      figure = create_figure(
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        dependencies: {
          vulcan: "sha256:abc",
        },
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)

      setup_workflow_mocks

      auth_header(:viewer)
      contents = {
        title: "Hercules Fighting the Nemean Lion",
        comment: "Better title",
        update_dependencies: true,
      }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        title: "Hercules Fighting the Nemean Lion",
        workflow_name: "test_workflow.cwl",
      )

      # Using the update_dependencies flag should update dependencies
      expect(json_body[:dependencies]).not_to eq({
        vulcan: "sha256:abc",
      })
      expect(Vulcan::WorkflowSnapshot.count).to eq(2)

      # Generally checks that metadata should update
      expect(Vulcan::WorkflowSnapshot.first.authors).not_to eq(
        Vulcan::WorkflowSnapshot.last.authors
      )

      # Make sure the workflow steps also update
      expect(Vulcan::WorkflowSnapshot.first.cwl_yaml).not_to eq(
        Vulcan::WorkflowSnapshot.last.cwl_yaml
      )
    end

    it "does not capture workflow updates during figure update, if dependencies not updated" do
      figure = create_figure(
        title: "Lion of Nemea",
        workflow_name: "test_workflow.cwl",
        dependencies: {
          vulcan: "sha256:abc",
        },
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)

      setup_workflow_mocks

      auth_header(:viewer)
      contents = {
        title: "Hercules Fighting the Nemean Lion",
        comment: "Better title",
      }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        title: "Hercules Fighting the Nemean Lion",
        workflow_name: "test_workflow.cwl",
        dependencies: {
          vulcan: "sha256:abc",
        },
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(2)
      expect(Vulcan::WorkflowSnapshot.first.cwl_yaml).to eq(
        Vulcan::WorkflowSnapshot.last.cwl_yaml
      )
    end

    it "guest updates an existing figure with non-public tags" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "reubens")
      auth_header(:guest)
      contents = { tags: ["private"], comment: "Make private" }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        tags: ["private"],
        workflow_name: "reubens",
      )
    end

    it "guest cannot update an existing figure with public tag" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")
      auth_header(:guest)
      contents = { tags: ["public"], comment: "Make public" }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(422)

      figure.refresh
      expect(figure.tags).to eq(nil)
    end
  end

  context "#delete" do
    it "deletes an existing figure" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(1)
      figure.refresh
      expect(figure.archived).to be_truthy
    end

    it "deletes an existing figure as guest" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:guest)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(1)
      figure.refresh
      expect(figure.archived).to be_truthy
    end

    it "throws exception for unknown figure" do
      figure = create_figure(title: "Lion of Nemea", workflow_name: "test_workflow.cwl")

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/99999999")

      expect(last_response.status).to eq(404)
      expect(Vulcan::Figure.count).to eq(1)
    end
  end
end
