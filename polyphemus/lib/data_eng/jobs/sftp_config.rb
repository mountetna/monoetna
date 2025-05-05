module WithSftpConfig
  def project_name
    config['project_name']
  end

  def workflow_config_id
    config['config_id']
  end

  def workflow_version
    config['version_number']
  end

  def magic_string
    config['config']['magic_string']
  end

  def name_regex
    Regexp.new("#{magic_string}(-|_)")
  end

  def file_regex
    Regexp.new("#{magic_string}(-|_).*")
  end

  def initial_start_scan_time
    runtime_config['config']['initial_start_scan_time']
  end

  def override_interval
    runtime_config['config']['override_interval']
  end

  def deposit_root_path
    # Mandatory params
    ::File.join('/', config['config']['deposit_root_path'] || '')
  end

  def deposit_host
    config["secrets"]["sftp_deposit_host"]
  end

  def deposit_user
    config["secrets"]["sftp_deposit_user"]
  end

  def deposit_password
    config["secrets"]["sftp_deposit_password"]
  end

  def override_root_path
    runtime_config['config']['override_root_path'].yield_self do |path|
      path == '' ? nil : path
    end
  end

  def raw_ingest_root_path
    config['config']['ingest_root_path']
  end

  def ingest_root_path
    override_root_path || raw_ingest_root_path
  end

  def ingest_host
    config["secrets"]["sftp_ingest_host"]
  end

  def ingest_user
    config["secrets"]["sftp_ingest_user"]
  end

  def ingest_password
    config["secrets"]["sftp_ingest_password"]
  end

  def metis_root_path
    config['config']['metis_root_path']
  end

  def bucket_name
    config['config']['bucket_name']
  end

  def notification_channel
    config["config"]["notification_channel"]
  end

  def notification_webhook_url
    config["secrets"]["notification_webhook_url"]
  end
end
