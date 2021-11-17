require "etna"
require_relative "../commands"

class Mvir1Waiver
  include WithEtnaClients
  include WithLogger

  def initialize(metis_client:, project_name: "mvir1")
    @metis_client = metis_client
    @project_name = project_name.to_s
    @release_bucket_name = bucket_name(:release_bucket)
    @restrict_bucket_name = bucket_name(:restrict_bucket)
  end

  def patient_name_search_string(patient_name)
    # This expects a Patient name like MVIR1-HS10,
    #   not an assay name, like MVIR1-HS10-D4BLD1-CYM2
    "#{patient_name}-%"
  end

  def pool_search_string(pool_name)
    # Slightly different regex for pools because we should be able
    #   to get the whole pool_name instead of a patient_name
    #   that we can't match to folder names.
    "#{pool_name}"
  end

  def diff_worker_for(release_folder, restricted_folder) end

  def empty_for(from_bucket, to_bucket, folder_name)
    Etna::Clients::Metis::WalkMetisDiffWorkflow.new(
      left_walker: Etna::Clients::Metis::WalkMetisWorkflow.new(
        metis_client: @metis_client,
        project_name: @project_name,
        bucket_name: from_bucket,
        root_dir: folder_name,
      ),
      right_walker: Etna::Clients::Metis::WalkMetisWorkflow.new(
        metis_client: @metis_client,
        project_name: @project_name,
        bucket_name: to_bucket,
        root_dir: folder_name,
      ),
    ).each do |diff|
      type, left, right = diff

      case type
      when :equal
        try_empty_file(left)
      when :left_unique
        try_copy(left, to_bucket)
        try_empty_file(left)
      when :right_unique
        # leave it be, maybe an artifact of concurrent behavior
      else
        try_remove(right)
        try_copy(left, to_bucket)
        try_empty_file(left)
      end
    end
  end

  # Try to atomically empty the contents of the file without it disappearing.
  def try_empty_file(file)
    return unless file.is_a?(Etna::Clients::Metis::File)

    uploader = Etna::Clients::Metis::MetisUploadWorkflow.new(
      metis_client: @metis_client,
      project_name: @project_name,
      bucket_name: file.bucket_name,
      metis_uid: '7af8e4be837744ce791c6239d00afefd',
    )

    uploader.do_upload(
      Etna::Clients::Metis::MetisUploadWorkflow::StreamingIOUpload.new(
        readable_io: StringIO.new
      ),
      file.file_path
    )
  end

  def try_copy(file, to_bucket)
    if file.is_a?(Etna::Clients::Metis::Folder)
      req = Etna::Clients::Metis::CreateFolderRequest.new(
        project_name: @project_name,
        bucket_name: to_bucket,
        folder_path: file.folder_path
      )
      @metis_client.create_folder(req) unless @metis_client.folder_exists?(req)
      return
    end

    downloader = Etna::Clients::Metis::MetisDownloadWorkflow.new(
      metis_client: @metis_client,
      project_name: @project_name,
      bucket_name: file.bucket_name
    )

    io = StringIO.new
    downloader.do_download(io, file)
    io.rewind

    uploader = Etna::Clients::Metis::MetisUploadWorkflow.new(
      metis_client: @metis_client,
      project_name: @project_name,
      bucket_name: to_bucket,
      metis_uid: '7af8e4be837744ce791c6239d00afefd',
    )

    uploader.do_upload(
      Etna::Clients::Metis::MetisUploadWorkflow::StreamingIOUpload.new(
        readable_io: io,
        size_hint: io.length,
      ),
      file.file_path
    )
  end

  def try_remove(file)
    if file.is_a?(Etna::Clients::Metis::File)
      @metis_client.delete_file(Etna::Clients::Metis::DeleteFileRequest.new(
        project_name: @project_name,
        bucket_name: file.bucket_name,
        file_path: file.file_path,
      ))
    end
  end

  def restrict_patient_data(patient_name, delete_on_metis)
    folders = find_folders_in_bucket(
      @release_bucket_name, patient_name_search_string(patient_name))

    if delete_on_metis
      # Completely move and remove source files.
      @metis_client.rename_folders(
        project_name: @project_name,
        source_bucket: @release_bucket_name,
        source_folders: folders,
        dest_bucket: @restrict_bucket_name,
      )
    else
      # 'copy' over files but leave empty marker files for anything missing.
      folders.each do |folder|
        empty_for(@release_bucket_name, @restrict_bucket_name, folder.folder_name)
      end
    end
  end

  def release_patient_data(patient_name)
    @metis_client.rename_folders(
      project_name: @project_name,
      source_bucket: @restrict_bucket_name,
      source_folders: find_folders_in_bucket(@restrict_bucket_name, patient_name_search_string(patient_name)),
      dest_bucket: @release_bucket_name,
    )
  end

  def restrict_pool_data(pool_name)
    @metis_client.rename_folders(
      project_name: @project_name,
      source_bucket: @release_bucket_name,
      source_folders: find_folders_in_bucket(@release_bucket_name, pool_search_string(pool_name)),
      dest_bucket: @restrict_bucket_name,
    )
  end

  def release_pool_data(pool_name)
    @metis_client.rename_folders(
      project_name: @project_name,
      source_bucket: @restrict_bucket_name,
      source_folders: find_folders_in_bucket(@restrict_bucket_name, pool_search_string(pool_name)),
      dest_bucket: @release_bucket_name,
    )
  end

  private

  def bucket_name(bucket_name)
    Polyphemus.instance.config(:metis)[bucket_name]
  end

  def find_folders_in_bucket(bucket_name, search_string)
    @metis_client.find(Etna::Clients::Metis::FindRequest.new(
      project_name: @project_name,
      bucket_name: bucket_name,
      params: [Etna::Clients::Metis::FindParam.new(
        attribute: "name",
        type: "folder",
        predicate: "=~",
        value: search_string,
      )],
    )).folders.all
  end
end
