require_relative '../lib/server/controllers/figure_controller'

describe Vulcan::Figure do
  before(:each) do
    clear_store
  end

  it 'returns a list of thumbnails' do
    store(
      'abee47d3ee8ba11e3fc2706d8d258e9e586465b4',
      'thumb.png',
      'thumbnail'
    )
    figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')

    expect(figure.thumbnails(storage: Vulcan::Storage.new)).to eq(['https://vulcan.test/api/labors/data/abee47d3ee8ba11e3fc2706d8d258e9e586465b4/thumb.png'])
  end

  it 'ignores unbuilt thumbnails' do
    figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')

    expect(figure.thumbnails(storage: Vulcan::Storage.new)).to eq([])
  end
end

describe FigureController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'returns a list of figures' do
      auth_header(:viewer)
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')
      get("/api/labors/figures")

      expect(last_response.status).to eq(200)
      expect(json_body[:figures].first[:title]).to eql(figure.title)
    end
  end

  context '#get' do
    it 'returns a figure by id' do
      auth_header(:viewer)
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')
      get("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(json_body[:title]).to eql(figure.title)
    end
  end

  context '#create' do
    it 'creates a new figure' do
      expect(Vulcan::WorkflowSnapshot.count).to eq(0)
      auth_header(:viewer)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        inputs: { a: 'b' }
      }
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors"
        )
      )

      expect(json_body.keys).to include(:dependencies)
      expect(Vulcan::WorkflowSnapshot.count).to eq(1)
    end

    it 'creates a new figure with tags' do
      auth_header(:viewer)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        inputs: { a: 'b' },
        tags: ['public']
      }
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors"
        )
      )

      expect(json_body.keys).to include(:dependencies, :workflow_snapshot)
      expect(json_body[:workflow_snapshot]).not_to eq(nil)
    end
  end

  context '#update' do
    it 'updates an existing figure' do
      figure = create_figure(
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        dependencies: {
          something: 'sha:abc'
        }
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)
      auth_header(:viewer)
      contents = { title: 'Hercules Fighting the Nemean Lion', comment: 'Better title'}
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        title: "Hercules Fighting the Nemean Lion",
        workflow_name: "test_workflow.cwl",

        # Updating maintains the previous figure dependencies
        dependencies: {
          something: 'sha:abc'
        }
      )

      # there are two figures
      expect(Vulcan::Figure.count).to eq(2)
      expect(Vulcan::WorkflowSnapshot.count).to eq(2)
      expect(Vulcan::WorkflowSnapshot.first.to_workflow_json).to eq(
        Vulcan::WorkflowSnapshot.last.to_workflow_json)
    end

    it 'throws exception for unknown figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')
      auth_header(:viewer)
      contents = { title: 'Hercules Fighting the Nemean Lion', comment: 'Better title' }
      post("/api/labors/figure/999999999/update", contents)

      expect(last_response.status).to eq(404)
    end

    it 'updates an existing figure with tags' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')
      auth_header(:viewer)
      contents = { tags: ['private'], comment: 'Make private' }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        tags: ['private'],
        workflow_name: "test_workflow.cwl"
      )
    end

    it 'can update dependencies' do
      figure = create_figure(
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        dependencies: {
          something: 'sha:abc'
        }
      )

      auth_header(:viewer)
      contents = {
        title: 'Hercules Fighting the Nemean Lion',
        comment: 'Better title',
        update_dependencies: true
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
        something: 'sha:abc'
      })

      # CI and local don't necessarily have production-tagged images, so
      #   the SHA's will just come back as empty strings.
      # But this is what should happen.
      # expect(json_body[:dependencies].values.all? { |d| d =~ /^sha256:.*/ }).to eq(true)
    end

    def setup_workflow_mocks
      allow(Etna::Cwl::Workflow).to receive(:metadata).and_return({
        "authors": ["a new author"],
        "projects": ["another project"]
      })
      allow_any_instance_of(Etna::Cwl::Workflow).to receive(:as_steps_json).and_return({
        "steps": [],
        "name": "test_workflow.cwl"
      })
      allow_any_instance_of(FigureController).to receive(:dependency_shas).and_return({
        "vulcan": "sha256:123"
      })
    end

    it 'captures any workflow updates during figure update, if dependencies updated' do
      figure = create_figure(
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        dependencies: {
          vulcan: 'sha256:abc'
        }
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)

      setup_workflow_mocks

      auth_header(:viewer)
      contents = {
        title: 'Hercules Fighting the Nemean Lion',
        comment: 'Better title',
        update_dependencies: true
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
        vulcan: 'sha256:abc'
      })
      expect(Vulcan::WorkflowSnapshot.count).to eq(2)

      # Generally checks that metadata should update
      expect(Vulcan::WorkflowSnapshot.first.authors).not_to eq(
        Vulcan::WorkflowSnapshot.last.authors)

      # Make sure the workflow steps also update
      expect(Vulcan::WorkflowSnapshot.first.steps).not_to eq(
        Vulcan::WorkflowSnapshot.last.steps)
    end

    it 'does not capture workflow updates during figure update, if dependencies not updated' do
      figure = create_figure(
        title: 'Lion of Nemea',
        workflow_name: 'test_workflow.cwl',
        dependencies: {
          vulcan: 'sha256:abc'
        }
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(1)

      setup_workflow_mocks

      auth_header(:viewer)
      contents = {
        title: 'Hercules Fighting the Nemean Lion',
        comment: 'Better title'
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
          vulcan: 'sha256:abc'
        }
      )

      expect(Vulcan::WorkflowSnapshot.count).to eq(2)
      expect(Vulcan::WorkflowSnapshot.first.to_workflow_json).to eq(
        Vulcan::WorkflowSnapshot.last.to_workflow_json)
    end
  end

  context '#delete' do
    it 'deletes an existing figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(1)
      figure.refresh
      expect(figure.archived).to be_truthy
    end

    it 'throws exception for unknown figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow.cwl')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/99999999")

      expect(last_response.status).to eq(404)
      expect(Vulcan::Figure.count).to eq(1)
    end
  end
end
