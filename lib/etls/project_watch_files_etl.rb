require_relative '../metis_file_in_watch_folder_etl'
require_relative './watch_etl_config'

class Polyphemus::ProjectWatchFilesEtl < Polyphemus::MetisFileInWatchFolderEtl
  attr_reader :config

  def initialize(config)
    @config = config
    super(bucket_watch_configs: config.bucket_configs, limit: 100)
  end

  def process_files(cursor, files, watch_type)
    bucket_name = cursor[:bucket_name]

    processors = config.file_processors(bucket_name, watch_type)
    return if processors.empty?
    processors.each { |p| p.process(cursor, files) }
  end
end
