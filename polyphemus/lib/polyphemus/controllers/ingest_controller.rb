require_relative "controller"

class IngestController < Polyphemus::Controller
  def list_dir
    require_params(:project_name, :folder_path, :ingest_host)
    require_valid_host

    files = Polyphemus::IngestFile.where(
      host: host,
      name: /^#{folder_path}/,
      ingested_at: nil,
      should_ingest: false,
    ).all

    success({ files: files.map { |r| r.to_hash } }.to_json, "application/json")
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end

  def queue
    require_params(:project_name, :folder_path, :ingest_host)
    require_valid_host

    files = Polyphemus::IngestFile.where(
      host: host,
      name: /^#{folder_path}/,
      ingested_at: nil,
      should_ingest: false,
    ).all do |file|
      file.update(should_ingest: true)
      file.refresh
    end

    success({ files: files.map { |r| r.to_hash } }.to_json, "application/json")
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end

  private

  def host
    @params[:ingest_host]
  end

  def folder_path
    path = @params[:folder_path]
    return path if "/" == path[-1]

    "#{path}/"
  end

  def require_valid_host
    raise Etna::BadRequest, "Invalid ingest host." unless configured_hosts.include?(host)
  end

  def configured_hosts
    Polyphemus.instance.config(:ingest)[:sftp].map { |c| c[:host] }
  end
end
