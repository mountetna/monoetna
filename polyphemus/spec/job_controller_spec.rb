describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it 'project editors cannot submit jobs' do
    post('/test/job', model_names: ["all"], redcap_tokens: ["123"])

    expect(last_response.status).to eq(401)
  end

  it 'model names must be an array' do
    auth_header(:administrator)
    post('/test/job', model_names: "all", redcap_tokens: ["123"])

    expect(last_response.status).to eq(422)
  end

  it 'redcap tokens must be an array' do
    auth_header(:administrator)
    post('/test/job', model_names: ["all"], redcap_tokens: "123")

    expect(last_response.status).to eq(422)
  end

  it 'project administrators can submit jobs' do
    stub_magma_models
    stub_magma_update_json
    stub_redcap_data

    auth_header(:administrator)
    post('/test/job', model_names: ["all"], redcap_tokens: ["123"])

    require 'pry'
    binding.pry
    expect(last_response.status).to eq(200)
    expect(last_response.body.include?(":model_one=>")).to eq(true)
    expect(last_response.body.include?(":model_two=>")).to eq(true)
  end
end

