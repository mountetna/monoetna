require_relative "./metis_folder_filtering_base_etl"
require_relative './watch_etl_config'

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
    @bucket_watch_configs.find { |b| b.bucket_name == cursor[:bucket_name] && b.project_name == cursor[:project_name] }.tap do |config|
      if config.nil?
        raise "Could not find bucket_config matching cursor, either reset the cursor or fix the config"
      end
    end
  end

  def watch_type_groups(cursor, folders)
    return {} if bucket_config(cursor).nil?

    {}.tap do |watch_type_groups|
      folders.each do |folder|
        bucket_config(cursor).find_matching_watches(folder.folder_path).each do |watch_config|
          f = watch_type_groups[watch_config.watch_type] ||= []
          f << folder
        end
      end
    end
  end

  def process(cursor, folders)
    # We'll need to filter out the folders based on folder_name here, using the supplied regexes.
    logger.info("Found #{folders.length} updated matching folders: #{folders.map { |f| f.folder_path }.join(",")}")

    target_folders = filter_target_folders(cursor, folders)
    logger.info("Found #{target_folders.length} folders that match the cursor regex: #{target_folders.map { |f| f.folder_path }.join(",")}")

    watch_type_groups(cursor, target_folders).each do |watch_type, watch_folders|
      process_watch_type_group(cursor, watch_type, watch_folders, folders)
    end

    logger.info("Done")
  end

  def process_watch_type_group(cursor, watch_type, watch_folders, all_batch_folders)
    current_watch_folders = self.current_watch_folders(cursor, watch_folders, watch_type)
    current_folder_metis_ids = current_watch_folders.map(&:metis_id)
    new_folders_to_watch, _ = partition_folders(watch_folders, current_folder_metis_ids)
    logger.info("Found #{new_folders_to_watch.length} new folders to watch as #{watch_type}: #{new_folders_to_watch.map { |f| f.folder_path }.join(",")}")


    process_folder_contents(cursor, watch_folders, watch_type)
    # Do not create the watch folders until we have successfully processed contents, else we risk 'missing' processing
    # a folder content due to an error between creating the records (which would prevent those folders from showing up
    # in the next new_folders_to_watch) and processing folder contents.
    # The guarantee: if we've created a watch folder record, it's because we also successfully processed its contents once.
    apply_records(cursor, watch_folders, watch_type)
    current_folder_metis_ids += new_folders_to_watch.map(&:id)
    # It's possible for dupes to exist in the case of an existing watch folder in the metis_ids
    # but it is newer than the cursor
    current_folder_metis_ids.uniq!
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
    Polyphemus.instance.db[:watch_folders].insert_conflict(
      target: [:bucket_name, :project_name, :metis_id, :watch_type],
      update: {
        folder_path: Sequel[:excluded][:folder_path],
        updated_at: Sequel[:excluded][:updated_at],
      }
    ).multi_insert(folders.map do |folder|
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

  def current_watch_folders(cursor, folders, watch_type)
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
    )

    find_request.add_param(Etna::Clients::Metis::FindParam.new(
      type: "file",
      attribute: "folder_id",
      predicate: "in",
      value: folders.map(&:id),
    ))

    metis_client.find(find_request).files.all
  end
end
