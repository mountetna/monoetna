require_relative "../metis_folder_etl"

class Polyphemus::ProjectPropagateFolderUpdatedAtEtl < Polyphemus::MetisFolderEtl
  def initialize(config, limit:)
    super(
      project_bucket_pairs: [config.project_name].product(config.buckets),
      limit: limit
    )
  end

  def process(cursor, folders)
    # For each folder that is updated, we'll touch it's child folders
    logger.info("Found #{folders.length} updated folders: #{folders.map { |f| f.folder_path }.join(",")}")
    update_child_folders(folders)
    logger.info("Done")
  end

  def update_child_folders(folders)
    # Only update if child folder's updated_at < parent folder's updated_at
    folders.each do |parent_folder|
      child_folders = children_folders(parent_folder)

      folders_to_touch = outdated_child_folders(parent_folder, child_folders)

      logger.info("Found #{folders_to_touch.length} child folders that will be touched: #{folders_to_touch.map { |f| f.folder_path }.join(",")}") unless folders_to_touch.empty?

      touch_folders(folders_to_touch)
    end
  end

  def children_folders(parent_folder)
    self.metis_client.list_folder_by_id(
      Etna::Clients::Metis::ListFolderByIdRequest.new(
        project_name: parent_folder.project_name,
        bucket_name: parent_folder.bucket_name,
        folder_id: parent_folder.id,
      )
    ).folders.all
  end

  def outdated_child_folders(parent_folder, child_folders)
    child_folders.select do |child|
      child.updated_at < parent_folder.updated_at
    end
  end

  def touch_folders(folders)
    folders.each do |folder|
      self.metis_client.touch_folder(
        Etna::Clients::Metis::TouchFolderRequest.new(
          project_name: folder.project_name,
          bucket_name: folder.bucket_name,
          folder_path: folder.folder_path,
        )
      )
    end
  end
end
