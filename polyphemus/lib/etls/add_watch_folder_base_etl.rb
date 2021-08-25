require_relative "../metis_folder_etl"

class Polyphemus::AddWatchFolderBaseEtl < Polyphemus::MetisFolderEtl
  def initialize(project_bucket_pairs:, model_name:, folder_name_globs: [], limit: 20)
    @model_name = model_name
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

    logger.info("Done")
  end

  def new_folders(cursor, folders)
    current_folder_paths = current_watch_folder_paths(cursor)

    folders.select do |folder|
      !(current_folder_paths.include?(folder.folder_path) || should_skip_folder?(folder))
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
      }
    end)
  end

  def should_skip_folder?(folder)
    # Subclasses can override
    false
  end

  def current_watch_folder_paths(cursor)
    Polyphemus::WatchFolder.where(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
    ).all.map do |folder|
      folder.folder_path
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
    require "pry"
    binding.pry
    # This could be overridden by subclasses
    self.metis_client.list_folder(Etna::Clients::Metis::ListFolderRequest.new(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      folder_path: folder.folder_path,
    )).files.all
  end
end
