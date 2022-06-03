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
    figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow')

    expect(figure.thumbnails(storage: Vulcan::Storage.new)).to eq(['https://vulcan.test/api/labors/data/abee47d3ee8ba11e3fc2706d8d258e9e586465b4/thumb.png'])
  end

  it 'ignores unbuilt thumbnails' do
    figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow')

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
      below_editor_roles.each do |role|
        auth_header(role)
        figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow')
        get("/api/labors/figures")

        expect(last_response.status).to eq(200)
        expect(json_body[:figures].first[:title]).to eql(figure.title)
      end
    end
  end

  context '#get' do
    it 'returns a figure by id' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'test_workflow')

      below_editor_roles.each do |role|
        auth_header(role)
        get("/api/labors/figure/#{figure.figure_id}")

        expect(last_response.status).to eq(200)
        expect(json_body[:title]).to eql(figure.title)
      end
    end
  end

  context '#create' do
    it 'creates a new figure' do
      below_editor_roles.each do |role|
        auth_header(role)
        contents = {
          title: 'Lion of Nemea',
          workflow_name: 'reubens',
          inputs: { a: 'b' }
        }
        post("/api/labors/figure/create", contents)

        expect(last_response.status).to eq(200)
        expect(json_body).to include(
          contents.merge(
            project_name: "labors"
          )
        )
      end
    end

    it 'creates a new public figure with tags' do
      auth_header(:viewer)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'reubens',
        inputs: { a: 'b' },
        tags: ['public']
      }

      expect(Vulcan::Figure.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors"
        )
      )
      expect(Vulcan::Figure.count).to eq(1)
    end

    it 'refuses to create a public figure with tags for guest' do
      auth_header(:guest)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'reubens',
        inputs: { a: 'b' },
        tags: ['public']
      }

      expect(Vulcan::Figure.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(422)

      expect(Vulcan::Figure.count).to eq(0)
    end

    it 'guest can create a figure with non-public tags' do
      auth_header(:guest)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'reubens',
        inputs: { a: 'b' },
        tags: ['mine']
      }

      expect(Vulcan::Figure.count).to eq(0)
      post("/api/labors/figure/create", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        contents.merge(
          project_name: "labors"
        )
      )

      expect(Vulcan::Figure.count).to eq(1)
    end
  end

  context '#update' do
    it 'updates an existing figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      original_title = 'Hercules Fighting the Nemean Lion'
      contents = { title: original_title, comment: 'Better title'}

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
          workflow_name: "reubens"
        )

        # there are more figures
        expect(Vulcan::Figure.count).to eq(2 + index)
      end
    end

    it 'throws exception for unknown figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:viewer)
      contents = { title: 'Hercules Fighting the Nemean Lion', comment: 'Better title' }
      post("/api/labors/figure/999999999/update", contents)

      expect(last_response.status).to eq(404)
    end

    it 'updates an existing figure with tags' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:viewer)
      contents = { tags: ['private'], comment: 'Make private' }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        tags: ['private'],
        workflow_name: "reubens"
      )
    end

    it 'guest updates an existing figure with non-public tags' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:guest)
      contents = { tags: ['private'], comment: 'Make private' }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        tags: ['private'],
        workflow_name: "reubens"
      )
    end

    it 'guest cannot update an existing figure with public tag' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:guest)
      contents = { tags: ['public'], comment: 'Make public' }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(422)
      
      figure.refresh
      expect(figure.tags).to eq(nil)
    end
  end

  context '#delete' do
    it 'deletes an existing figure as viewer' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(1)
      figure.refresh
      expect(figure.archived).to be_truthy
    end

    it 'deletes an existing figure as guest' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:guest)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(1)
      figure.refresh
      expect(figure.archived).to be_truthy
    end

    it 'throws exception for unknown figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/99999999")

      expect(last_response.status).to eq(404)
      expect(Vulcan::Figure.count).to eq(1)
    end
  end
end
