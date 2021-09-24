require_relative "./metis_folder_filtering_base_etl"

class Polyphemus::AddWatchFolderBaseEtl < Polyphemus::MetisFolderFilteringBaseEtl
  def initialize(project_bucket_pairs:, model_name:, folder_path_regexes: {}, limit: 20, watch_type:, linker:)
    @linker = linker
    @model_name = model_name
    @watch_type = watch_type
    super(
      project_bucket_pairs: project_bucket_pairs,
      limit: limit,
      folder_path_regexes: folder_path_regexes,
    )
  end

  def process(cursor, folders)
    # We'll need to filter out the folders based on folder_name here, using the
    #   supplied regexes
    logger.info("Found #{folders.length} updated folders: #{folders.map { |f| f.folder_path }.join(",")}")

    target_folders = filter_target_folders(cursor, folders)

    logger.info("Found #{target_folders.length} folders that match the cursor regex: #{target_folders.map { |f| f.folder_path }.join(",")}")

    new_folders_to_watch = new_folders(cursor, target_folders)

    logger.info("Found #{new_folders_to_watch.length} new folders to watch: #{new_folders_to_watch.map { |f| f.folder_path }.join(",")}")

    create_records(cursor, new_folders_to_watch)
    process_folder_contents(cursor, new_folders_to_watch)
    remove_changed_records(cursor, folders)

    logger.info("Done")
  end

  def new_folders(cursor, folders)
    current_folder_metis_ids = current_watch_folder_metis_ids(cursor, folders)

    folders.select do |folder|
      !(current_folder_metis_ids.include?(folder.id) || should_skip_folder?(folder))
    end
  end

  def remove_changed_records(cursor, folders)
    current_folder_metis_ids = current_watch_folder_metis_ids(cursor, folders)

    folders_to_remove = folders.select do |folder|
      current_folder_metis_ids.include?(folder.id) && !folder_path_satisfies_regex?(cursor, folder)
    end

    logger.info("Found #{folders_to_remove.length} changed folders. Removing as watch folders: #{folders_to_remove.map { |f| f.folder_path }.join(",")}")

    Polyphemus::WatchFolder.where(
      metis_id: folders_to_remove.map { |f| f.id },
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      watch_type: @watch_type,
    ).all do |f|
      f.delete
    end
  end

  def create_records(cursor, folders)
    Polyphemus::WatchFolder.multi_insert(folders.map do |folder|
      {
        created_at: DateTime.now,
        updated_at: DateTime.now,
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        folder_path: folder.folder_path,
        watch_type: @watch_type,
        metis_id: folder.id,
      }
    end)
  end

  def current_watch_folder_metis_ids(cursor, folders)
    Polyphemus::WatchFolder.where(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      metis_id: folders.map do |folder|
        folder.id
      end,
      watch_type: @watch_type,
    ).all.map do |folder|
      folder.metis_id
    end
  end

  def process_folder_contents(cursor, folders)
    # This could be overridden by subclasses
    folders.each do |folder|
      files = folder_files(cursor, folder)
      @linker.link(
        model_name: @model_name,
        files: files,
      ) if @linker
    end
  end

  def folder_files(cursor, folder)
    # This could be overridden by subclasses
    self.metis_client.list_folder_by_id(Etna::Clients::Metis::ListFolderByIdRequest.new(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      folder_id: folder.id,
    )).files.all
  end
end
