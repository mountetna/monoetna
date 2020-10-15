describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it 'configuration fails for non-superusers' do
    auth_header(:viewer)
    get('/configuration')

    expect(last_response.status).to eq(401)
  end

  it 'superusers can fetch app configuration' do
    auth_header(:viewer)
    get('/configuration')

    expect(last_response.status).to eq(200)
    expect(json_body.keys).to eq('magma', 'metis', 'janus', 'timur', 'polyphemus')
  end

  it 'shows polyphemus is available for users' do
    auth_header(:viewer)
    get('/')

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('Polyphemus is available.')
  end
end

