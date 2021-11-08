module Polyphemus::WithIngestConfig
  def cat_config
    @cat_config ||= sftp_configs.find { |c| c[:alias] == "cat" }
  end

  def c4_config
    @c4_config ||= ssh_configs.find { |c| c[:alias] == "c4" }
  end

  def sftp_configs
    ingest_configs[:sftp]
  end

  def ssh_configs
    ingest_configs[:ssh]
  end

  def ingest_configs
    Polyphemus.instance.config(:ingest)
  end
end
