require_relative './add_watch_folder_base_etl'

# For now, this class is instantiated within ruby and managed in a static way.
# In the future, this object could be constructed by the executor based on project
# configuration, or pulled from several sources, or like whatever.  Doesn't matter.
# Point is, this is a good place for us to start modeling what project configuration
# looks like.
class Polyphemus::ProjectWatchFoldersConfig
  attr_reader :project_name

  def initialize(project_name:)
    @project_name = project_name
    @buckets = {}
    @processors = {}
  end

  def bucket(bucket_name)
    @buckets[bucket_name] ||= Polyphemus::BucketWatchFoldersConfig.new(bucket_name: bucket_name, project_name: project_name)
  end

  def buckets
    @buckets.keys
  end

  def bucket_configs
    @buckets.values
  end

  # processors are objects receiving process(cursor, files)
  def process_watch_type_with(watch_type_config, *processors)
    (@processors[watch_type_config.truple] ||= []).push(*processors)
    watch_type_config
  end

  def processors(bucket_name, watch_type)
    watch_type = bucket(bucket_name).watcher(watch_type)
    @processors[watch_type.truple] ||= []
  end
end


class Polyphemus::ProjectWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  attr_reader :config

  def initialize(config)
    @config = config
    super(bucket_watch_configs: config.bucket_configs)
  end

  def process_folder_contents(cursor, folders, watch_type)
    bucket_name = cursor[:bucket_name]

    processors = config.processors(bucket_name, watch_type)
    return if processor.empty?
    files = folders_files(cursor, folders)
    processors.each { |p| p.process(cursor, files) }
  end
end
