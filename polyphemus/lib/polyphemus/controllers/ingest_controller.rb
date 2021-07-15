require_relative "controller"

class IngestController < Polyphemus::Controller
  def list_dir
    require_params(:project_name, :folder_path, :ingest_host)
    require_valid_host

    results = Polyphemus::IngestFile.where(
      host: @params[:ingest_host],
      name: /^#{@params[:folder_path]}/,
      ingested_at: nil,
    ).all

    success({ results: results.map { |r| r.to_hash } }.to_json, "application/json")
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end

  def update
    require_params(:project_name, :folder_path, :ingest_host, :updates)
    require_valid_host

    success({ success: true }.to_json, "application/json")
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end

  private

  def require_valid_host
    raise Etna::BadRequest, "Invalid ingest host." unless configured_hosts.include?(@params[:ingest_host])
  end

  def configured_hosts
    Polyphemus.instance.config(:ingest)[:sftp].map { |c| c[:host] }
  end
end
