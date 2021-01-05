describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it 'project editors cannot submit jobs' do
    post('/test/job', model_names: ["all"], redcap_tokens: ["123"])

    expect(last_response.status).to eq(401)
  end

  it 'model names must be an array or "all"' do
    auth_header(:administrator)
    post('/test/job', model_names: "blah", redcap_tokens: ["123"])

    expect(last_response.status).to eq(422)
  end

  it 'redcap tokens must be an array' do
    auth_header(:administrator)
    post('/test/job', model_names: "all", redcap_tokens: "123")

    expect(last_response.status).to eq(422)
  end

  it 'project administrators can submit jobs' do
    # Not a great test ... can't figure out how to test or mock for
    #   a process spun out in a different Thread.
    stub_magma_models
    stub_magma_update_json
    stub_redcap_data

    auth_header(:administrator)
    post('/test/job', model_names: "all", redcap_tokens: ["123"])

    expect(last_response.status).to eq(200)

    expect(json_body[:record].keys.first).to eq(:model_one)
    expect(json_body[:record].keys.last).to eq(:model_two)
  end
end

