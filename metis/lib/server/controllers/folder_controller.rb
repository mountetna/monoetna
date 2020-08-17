require_relative '../../folder_rename_revision'

class FolderController < Metis::Controller
  def list
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    files = Metis::File.where(
      project_name: @params[:project_name],
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |file|
      file.to_hash(@request)
    end

    folders = Metis::Folder.where(
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |fold|
      fold.to_hash
    end

    success_json(files: files, folders: folders)
  end

  def list_all_folders
    bucket = require_bucket

    # We start with the root folders, and iterate
    #   from there.
    parent_folder_ids = nil
    folders = []
    have_leaves = true

    loop do
      child_folders = Metis::Folder.where(
        bucket: bucket,
        folder_id: parent_folder_ids
      ).all

      break if child_folders.length == 0

      parent_folder_ids = child_folders.map { |fold| fold[:id] }
      folders += child_folders
    end

    success_json(folders: folder_hashes_with_calculated_paths(folders))
  end

  def create
    bucket = require_bucket
    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:folder_path])

    folders = mkdir_p(bucket, @params[:folder_path], @params[:project_name], Metis::File.author(@user))
    success_json(folders: [ folders.first.to_hash ])
  end

  def remove
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, 'Folder is read-only' if folder.read_only?

    raise Etna::Forbidden, 'Parent folder is read-only' if folder.folder&.read_only?

    raise Etna::BadRequest, 'Folder is not empty' unless folder.can_remove?

    response = { folders: [ folder.to_hash ] }

    # actually remove the folder
    folder.remove!

    return success_json(response)
  end

  def protect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    folder.protect!

    success_json(folders: [ folder.to_hash ])
  end

  def unprotect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::BadRequest, 'Folder is not protected' unless folder.read_only?

    folder.unprotect!

    success_json(folders: [ folder.to_hash ])
  end

  def rename
    require_param(:new_folder_path)

    source_bucket = require_bucket

    dest_bucket = source_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    revision = Metis::FolderRenameRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:folder_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        dest_bucket.name,
        @params[:new_folder_path]
      ),
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    # we need the dest parent folder here, so we call File#path_parts
    source_folder_path, _ = Metis::File.path_parts(@params[:folder_path])
    dest_folder_path, _ = Metis::File.path_parts(@params[:new_folder_path])

    revision.set_folder(
      revision.source,
      [require_folder(source_bucket, source_folder_path)].compact)
    revision.set_folder(
      revision.dest,
      [require_folder(dest_bucket, dest_folder_path)].compact) if dest_folder_path

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    return success_json(folders: [ revision.revise!.to_hash ])
  end

  protected

  def mkdir_p(bucket, folder_path, project_name, author)
    existing_folders = Metis::Folder.from_path(bucket, folder_path, allow_partial_match: true)
    folder_names = folder_path.split(%r!/!)

    folder_names.inject([]) do |parents, folder_name|
      existing = existing_folders.shift
      next (parents << existing) unless existing.nil?


      if Metis::File.exists?(folder_name, bucket, parents.last)
        raise Etna::BadRequest, "Cannot overwrite existing file"
      end

      begin
        parents << Metis::Folder.create(
          folder_name: folder_name,
          bucket_id: bucket&.id,
          folder_id: parents.last&.id,
          project_name: project_name,
          author: author,
      )
      rescue  Sequel::UniqueConstraintViolation => e
        ## Can occur if two simult requests get to the create line after both reading no existing_folders.
        ## Because of the uniq index constraint, this will occur for one of the requests, while the other should succeed.
        ## In that case, fall back to querying the other transaction's entry.
        ## Note: find_or_create does not fix this, it still does not handle the unique constraint and simply
        ## queries or creates, which is not good enough for READ COMMITTED isolation where a read might not see a yet
        ## committed create.
        Metis.instance.logger.log_error(e)
        parents << Metis::Folder.find(bucket_id: bucket&.id, folder_id: parents.last&.id, folder_name: folder_name)
      end
    end
  end

  def folder_hashes_with_calculated_paths(folders)
    # Calculate the folder_path, instead of
    #   doing it in the database

    # To reduce database calls for each folder, we'll
    #   store a Hash of folder_id: path, and construct
    #   the folder_path attribute as we iterate through the folders.
    path_map = {}

    folders.map do |fold|
      folder_id_sym = fold.id.to_s.to_sym

      path_map[folder_id_sym] = fold.folder_id ?
        "#{path_map[fold.folder_id.to_s.to_sym]}/#{fold.folder_name}" :
        fold.folder_name
      folder_hash = fold.to_hash(false)
      folder_hash[:folder_path] = path_map[folder_id_sym]
      folder_hash
    end
  end
end
