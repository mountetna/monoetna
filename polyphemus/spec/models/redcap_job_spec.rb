describe Polyphemus::RedcapJob do
  def create_job(request_params: {job_type: 'redcap', project_name: PROJECT}, request_env: {}, response: {}, user: {})
    Polyphemus::RedcapJob.new(
      request_params: request_params,
      token: user[:token])
  end

  def payload_stamp
    {
      job_type: 'redcap',
      project_name: PROJECT,
      mode: 'default'
    }
  end

  context 'raises exception' do
    it 'if missing required job parameters' do
      redcap_job = create_job(
        request_params: payload_stamp.update(
          redcap_tokens: REDCAP_TOKEN
        )
      )

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "request_params missing required param(s): model_names")

      redcap_job = create_job(request_params: payload_stamp.update(model_names: "all"))

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "request_params missing required param(s): redcap_tokens")
    end

    it 'if model_names param is invalid' do
      redcap_job = create_job(request_params: payload_stamp.update(redcap_tokens: REDCAP_TOKEN, model_names: "Foo"))

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "model_names must be \"all\" or a comma-separated list of model names.")

      redcap_job = create_job(request_params: payload_stamp.update(redcap_tokens: REDCAP_TOKEN, model_names: {one: "two", three: "four"}))

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "model_names must be \"all\" or a comma-separated list of model names.")
    end

    it 'if redcap_tokens param is invalid' do
      redcap_job = create_job(request_params: payload_stamp.update(redcap_tokens: [REDCAP_TOKEN], model_names: "all"))

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "redcap_tokens must be a comma-separated list of tokens.")

      redcap_job = create_job(request_params: payload_stamp.update(redcap_tokens: {project_one: "123", project_two: "345"}, model_names: "all"))

      expect {
        redcap_job.validate
      }.to raise_error(Polyphemus::JobError, "redcap_tokens must be a comma-separated list of tokens.")
    end
  end
end
