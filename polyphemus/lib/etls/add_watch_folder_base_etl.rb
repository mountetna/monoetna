require_relative "./metis_folder_filtering_base_etl"

class Polyphemus::BucketWatchFoldersConfig
  attr_reader :bucket_name, :project_name

  def initialize(bucket_name:, project_name:)
    @bucket_name = bucket_name
    @project_name = project_name
    @watches = {}
  end

  def cursor_pair
    [project_name, bucket_name]
  end

  def watcher(watch_type)
    @watches[watch_type] ||= Polyphemus::WatchTypeConfig.new(
      watch_type: watch_type,
      bucket_name: bucket_name,
      project_name: project_name
    )
  end

  def find_matching_watches(path)
    @watches.values.select do |config|
      path =~ config.sum_watch_regex
    end
  end

  def sum_bucket_regex
    Regexp.union(@watches.values.map(&:sum_watch_regex))
  end
end

class Polyphemus::WatchTypeConfig
  attr_reader :watch_type, :project_name, :bucket_name, :regexes

  def initialize(project_name:, bucket_name:, watch_type:)
    @watch_type = watch_type
    @project_name = project_name
    @bucket_name = bucket_name
    @regexes = []
  end

  def truple
    [project_name, bucket_name, watch_type]
  end

  def watch(regex)
    @regexes << regex unless @regexes.include?(regex)
    # Reset so the next sum call recalculates
    @sum_watch_regex = nil
    self
  end

  def sum_watch_regex
    @sum_watch_regex ||= Regexp.union(regexes)
  end
end

class Polyphemus::AddWatchFolderBaseEtl < Polyphemus::MetisFolderFilteringBaseEtl
  def initialize(bucket_watch_configs: [], limit: 20)
    pairs = bucket_watch_configs.map(&:cursor_pair)
    raise "bucket_watch_configs have non unique project/bucket combinations" unless pairs.uniq.length == pairs.length

    @bucket_watch_configs = bucket_watch_configs
    super(
      project_bucket_pairs: pairs,
      limit: limit,
      folder_path_regexes: bucket_watch_configs.map do |config|
        [
          project_bucket_symbol(
            project_name: config.project_name,
            bucket_name: config.bucket_name,
          ),
          config.sum_bucket_regex,
        ]
      end.to_h
    )
  end

  def bucket_config(cursor)
    @bucket_watch_configs.find { |b| b.bucket_name == cursor[:bucket_name] && b.project_name == cursor[:project_name] }
  end

  def watch_type_groups(cursor, folders)
    return {} if bucket_config.nil?

    {}.tap do |watch_type_groups|
      folders.each do |folder|
        bucket_config.find_matching_watches(folder).each do |watch_config|
          folders = watch_type_groups[watch_config.watch_type] ||= []
          folders << folder
        end
      end
    end
  end

  def process(cursor, folders)
    # We'll need to filter out the folders based on folder_name here, using the
    #   supplied regexes
    logger.info("Found #{folders.length} updated matching folders: #{folders.map { |f| f.folder_path }.join(",")}")

    target_folders = filter_target_folders(cursor, folders)
    logger.info("Found #{target_folders.length} folders that match the cursor regex: #{target_folders.map { |f| f.folder_path }.join(",")}")

    watch_type_groups(cursor, target_folders).each do |watch_type, watch_folders|
      process_watch_type_group(cursor, watch_type, watch_folders, folders)
    end

    logger.info("Done")
  end

  def process_watch_type_group(cursor, watch_type, watch_folders, all_batch_folders)
    current_folder_metis_ids = current_watch_folder_metis_ids(cursor, watch_folders, watch_type)
    new_folders_to_watch, _ = partition_folders(watch_folders, current_folder_metis_ids)
    logger.info("Found #{new_folders_to_watch.length} new folders to watch as #{watch_type}: #{new_folders_to_watch.map { |f| f.folder_path }.join(",")}")


    process_folder_contents(cursor, new_folders_to_watch, watch_type)
    # Do not create the watch folders until we have successfully processed contents, else we risk 'missing' processing
    # a folder content due to an error between creating the records (which would prevent those folders from showing up
    # in the next new_folders_to_watch) and processing folder contents.
    # The guarantee: if we've created a watch folder record, it's because we also successfully processed its contents once.
    apply_records(cursor, watch_folders, watch_type)
    current_folder_metis_ids += new_folders_to_watch.map(&:id)
    remove_changed_records(cursor, all_batch_folders, current_folder_metis_ids, watch_type)
  end

  def partition_folders(folders, current_folder_metis_ids)
    folders.partition do |folder|
      !current_folder_metis_ids.include?(folder.id)
    end
  end

  def remove_changed_records(cursor, folders, current_folder_metis_ids, watch_type)
    folders_to_remove = folders.select do |folder|
      current_folder_metis_ids.include?(folder.id) && !folder_path_satisfies_regex?(cursor, folder)
    end

    logger.info("Found #{folders_to_remove.length} changed folders. Removing as watch folders: #{folders_to_remove.map { |f| f.folder_path }.join(",")}")

    Polyphemus::WatchFolder.where(
      metis_id: folders_to_remove.map { |f| f.id },
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      watch_type: watch_type,
    ).all do |f|
      f.delete
    end
  end

  def apply_records(cursor, folders, watch_type)
    Polyphemus::WatchFolder.on_duplicate_key_update(:folder_path, :updated_at).multi_insert(folders.map do |folder|
      {
        created_at: DateTime.now,
        updated_at: DateTime.now,
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        folder_path: folder.folder_path,
        watch_type: watch_type,
        metis_id: folder.id,
      }
    end)
  end

  def current_watch_folder_metis_ids(cursor, folders, watch_type)
    Polyphemus::WatchFolder.where(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      metis_id: folders.map do |folder|
        folder.id
      end,
      watch_type: watch_type,
    ).all.map do |folder|
      folder.metis_id
    end
  end

  def process_folder_contents(cursor, folders, watch_type)
    # This could be overridden by subclasses to provide default behavior on new
    # folders that match a given watch.
  end

  def folders_files(cursor, folders)
    find_request = Etna::Clients::Metis::FindRequest.new(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      hide_paths: true,
    )

    find_request.add_param(Etna::Clients::Metis::FindParam.new(
      type: "file",
      attribute: "folder_id",
      predicate: "in",
      value: folders.map(&:id),
    ))

    metis_client.find(find_request).files.all.map(&:with_containing_folder)
  end
end
