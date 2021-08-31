require_relative "../metis_folder_etl"

class Polyphemus::AddWatchFolderBaseEtl < Polyphemus::MetisFolderEtl
  def initialize(project_bucket_pairs:, model_name:, folder_name_globs: [], limit: 20, watch_type:)
    @model_name = model_name
    @watch_type = watch_type
    super(
      project_bucket_pairs: project_bucket_pairs,
      folder_name_globs: folder_name_globs,
      limit: limit,
    )
  end

  def process(cursor, folders)
    # All of these folders should contain processed files that
    #   we want to watch and link in when / if they change.
    logger.info("Found #{folders.length} updated folders: #{folders.map { |f| f.folder_path }.join(",")}")

    new_folders_to_watch = new_folders(cursor, folders)

    logger.info("Found #{new_folders_to_watch.length} folders to watch: #{new_folders_to_watch.map { |f| f.folder_path }.join(",")}")

    create_records(cursor, new_folders_to_watch)
    process_folder_contents(cursor, new_folders_to_watch)
    remove_changed_records(cursor, folders)

    logger.info("Done")
  end

  def new_folders(cursor, folders)
    current_folder_ids = current_watch_folder_ids(cursor, folders)

    folders.select do |folder|
      !(current_folder_ids.include?(folder.id) || should_skip_folder?(folder))
    end
  end

  def remove_changed_records(cursor, folders)
    folders_to_remove = folders.select do |folder|
      !folder_path_satisfies_globs?(folder)
    end

    logger.info("Found #{folders_to_remove} changed folders. Removing as watch folders: #{folders_to_remove.map { |f| f.id }.join(",")}")

    Polyphemus::WatchFolder.where(
      folder_id: folders_to_remove.map { |f| f.id },
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
        folder_id: folder.id,
      }
    end)
  end

  def should_skip_folder?(folder)
    # Subclasses can override
    false
  end

  def folder_path_satisfies_globs?(folder)
    @folder_name_globs.all? do |glob|
      folder_subdirectories(folder.folder_path).any? do |path|
        ::File.fnmatch?(glob, path)
      end
    end
  end

  def folder_subdirectories(folder_path)
    # The Metis find API matches globs on subdirectories,
    #   so we'll have to test against those here.
    # For example, a "glob" of foo/* will return matches
    #   like:
    #
    #     - foo/sub-foo
    #     - parent-foo/foo/sub-foo
    parts = folder_path.split("/")
    [].tap do |results|
      (1..parts.length).each do
        results << parts.join("/")
        parts.shift
      end
    end
  end

  def current_watch_folder_ids(cursor, folders)
    Polyphemus::WatchFolder.where(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      folder_id: folders.map do |folder|
        folder.id
      end,
      watch_type: @watch_type,
    ).all.map do |folder|
      folder.folder_id
    end
  end

  def process_folder_contents(cursor, folders)
    # This could be overridden by subclasses
    folders.each do |folder|
      files = folder_files(cursor, folder)
      @linker.link(
        project_name: cursor[:project_name],
        model_name: @model_name,
        files: files,
      ) if @linker
    end
  end

  def folder_files(cursor, folder)
    # This could be overridden by subclasses
    self.metis_client.list_folder(Etna::Clients::Metis::ListFolderRequest.new(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      folder_path: folder.folder_path,
    )).files.all
  end
end
