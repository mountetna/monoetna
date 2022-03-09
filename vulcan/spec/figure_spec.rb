require_relative '../lib/server/controllers/figure_controller'

describe FigureController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'returns a list of figures' do
      auth_header(:viewer)
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      get("/api/labors/figures")

      expect(last_response.status).to eq(200)
      expect(json_body[:figures].first[:title]).to eql(figure.title)
    end
  end

  context '#get' do
    it 'returns a figure by id' do
      auth_header(:viewer)
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      get("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(json_body[:title]).to eql(figure.title)
    end
  end
  context '#create' do
    it 'creates a new figure' do
      auth_header(:viewer)
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

    it 'creates a new figure with tags' do
      auth_header(:viewer)
      contents = {
        title: 'Lion of Nemea',
        workflow_name: 'reubens',
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
    end
  end

  context '#update' do
    it 'updates an existing figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:viewer)
      contents = { title: 'Hercules Fighting the Nemean Lion' }
      post("/api/labors/figure/#{figure.figure_id}/update", contents)

      expect(last_response.status).to eq(200)
      expect(json_body).to include(
        figure_id: 1,
        inputs: {},
        project_name: "labors",
        title: "Hercules Fighting the Nemean Lion",
        workflow_name: "reubens"
      )
    end

    it 'throws exception for unknown figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:viewer)
      contents = { title: 'Hercules Fighting the Nemean Lion' }
      post("/api/labors/figure/999999999/update", contents)

      expect(last_response.status).to eq(404)
    end

    it 'updates an existing figure with tags' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')
      auth_header(:viewer)
      contents = { tags: ['private'] }
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
  end

  context '#delete' do
    it 'deletes an existing figure' do
      figure = create_figure(title: 'Lion of Nemea', workflow_name: 'reubens')

      expect(Vulcan::Figure.count).to eq(1)
      auth_header(:viewer)
      delete("/api/labors/figure/#{figure.figure_id}")

      expect(last_response.status).to eq(200)
      expect(Vulcan::Figure.count).to eq(0)
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
