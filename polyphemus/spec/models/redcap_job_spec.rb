# require_relative '../lib/folder_path_calculator'

describe Polyphemus::RedcapJob do
  def create_job(request_params: {job_type: 'redcap', project_name: PROJECT}, request_env: {}, response: {}, user: {})
    Polyphemus::RedcapJob.new(
      request_params: request_params,
      request_env: request_env,
      response: response,
      user: user)
  end

  def payload_stamp
    {
      job_type: 'redcap',
      project_name: PROJECT
    }
  end

  context 'raises exception' do
    it 'if missing required job parameters' do
      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        redcap_tokens: ["123"]
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "job_params missing required param(s): model_names")

      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        model_names: "all"
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "job_params missing required param(s): redcap_tokens")
    end

    it 'if model_names param is invalid' do
      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        redcap_tokens: ["123"],
        model_names: "foo"
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "model_names must be \"all\" or an array of model names.")

      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        redcap_tokens: ["123"],
        model_names: {one: "two", three: "four"}
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "model_names must be \"all\" or an array of model names.")
    end

    it 'if redcap_tokens param is invalid' do
      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        redcap_tokens: "123",
        model_names: "all"
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "redcap_tokens must be an array of tokens.")

      redcap_job = create_job(request_params: payload_stamp.update(job_params: {
        redcap_tokens: {project_one: "123", project_two: "345"},
        model_names: "all"
      }))

      expect {
        redcap_job.validate
      }.to raise_error(Etna::BadRequest, "redcap_tokens must be an array of tokens.")
    end
  end
end