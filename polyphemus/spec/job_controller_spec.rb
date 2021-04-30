describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before do
    copy_redcap_project
  end

  it 'project editors cannot submit jobs' do
    post('/test/job', job_type: "redcap", job_params: {
      model_names: ["all"], redcap_tokens: ["123"]
    })

    expect(last_response.status).to eq(401)
  end

  it 'model names must be an array or "all"' do
    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: "blah", redcap_tokens: ["123"]
    })

    expect(last_response.status).to eq(422)
  end

  it 'redcap tokens must be an array' do
    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: "all", redcap_tokens: "123"
    })

    expect(last_response.status).to eq(422)
  end

  it 'throws exception for unknown job type' do
    auth_header(:administrator)
    post('/test/job', job_type: "unsupported-job-type", job_params: {
      something: "important"
    })

    expect(last_response.status).to eq(422)
  end

  it 'throws exception for missing job type parameter' do
    auth_header(:administrator)
    post('/test/job', job_params: {
      something: "important"
    })

    expect(last_response.status).to eq(422)
  end

  it 'throws exception for missing job params parameter' do
    auth_header(:administrator)
    post('/test/job', job_type: "redcap")

    expect(last_response.status).to eq(422)
  end

  it 'throws exception for bad mode parameter' do
    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: "all", redcap_tokens: ["123"], mode: "rocks"
    })

    expect(last_response.status).to eq(422)
  end

  it 'project administrators can submit jobs' do
    # Not a great test ... can't figure out how to test or mock for
    #   a process spun out in a different Thread.
    stub_magma_models
    stub_magma_update_json
    stub_redcap_data

    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: "all", redcap_tokens: ["123"]
    })

    expect(last_response.status).to eq(200)

    expect(json_body[:results].keys).to match_array([:citation, :model_one, :model_two, :stats])

    # Updates all records found in REDCap, by default
    expect(json_body[:results][:model_two].keys).to match_array([:"123", :"321", :abc])
  end

  it 'only updates existing magma records in existing mode' do
    # Not a great test ... can't figure out how to test or mock for
    #   a process spun out in a different Thread.
    stub_magma_models
    stub_magma_update_json
    stub_redcap_data

    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: ["model_two"],
      redcap_tokens: ["123"],
      mode: "existing"
    })

    expect(last_response.status).to eq(200)

    expect(json_body[:results].keys).to eq([:model_two])
    expect(json_body[:results][:model_two].keys).to eq([:"123"])
  end

  it 'blanks "removed" magma records in strict mode' do
    # Not a great test ... can't figure out how to test or mock for
    #   a process spun out in a different Thread.
    stub_magma_models
    stub_magma_update_json
    stub_redcap_data

    auth_header(:administrator)
    post('/test/job', job_type: "redcap", job_params: {
      model_names: ["stats"],
      redcap_tokens: ["123"],
      mode: "strict"
    })

    expect(last_response.status).to eq(200)

    id_abc = temp_id(json_body[:results], "abc")

    expect(json_body[:results].keys).to match_array([:stats, :model_two])
    expect(json_body[:results][:model_two].keys).to match_array([:"123", :abc])
    expect(json_body[:results][:stats].keys).to match_array([id_abc])

    expect(json_body[:results][:model_two][:abc][:stats]).not_to eq([])
    expect(json_body[:results][:model_two][:"123"][:stats]).to eq([])
  end
end
