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
  end
end
