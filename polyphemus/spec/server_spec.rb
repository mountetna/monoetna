describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it 'configuration open to anyone' do
    get('/configuration')

    expect(last_response.status).to eq(200)
    expect(json_body[:test].keys.sort).to eq([:magma, :metis, :janus, :timur, :polyphemus, :auth_redirect, :docker].sort)
  end

  it 'superusers can fetch app configuration' do
    auth_header(:superuser)
    get('/configuration')

    expect(last_response.status).to eq(200)
    expect(json_body[:test].keys.sort).to eq([:magma, :metis, :janus, :timur, :polyphemus, :auth_redirect, :docker].sort)
  end

  it 'returns the polyphemus client' do
    auth_header(:viewer)
    get('/')

    expect(last_response.status).to eq(200)
    expect(last_response.body).to match('var CONFIG = ')
  end

  it 'returns a list of etl configs for a project' do
    auth_header(:editor)
    get('/api/labors/etl/configs')

    expect(last_response.status).to eq(200)
    expect(json_body.first.keys).to eq([:project_name, :etl, :name, :updated_at, :created_at, :config])
  end

  it 'returns a list of etl jobs for a project' do
    auth_header(:editor)
    get('/api/labors/etl/jobs')

    expect(last_response.status).to eq(200)
    expect(json_body.first.keys).to eq([:name, :schema])
  end
end

