class Metis
  class Folder < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :bucket

    many_to_one :folder

    one_to_many :folders

    one_to_many :files

    class CopyError < Exception
    end

    def self.mkdir_p(bucket, folder_path, project_name, author)
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

    def self.from_path(bucket, folder_path, allow_partial_match: false)
      return [] unless folder_path && !folder_path.empty? && bucket

      folder_names = folder_path.split(%r!/!)

      db = Metis.instance.db

      query = %Q(
        WITH RECURSIVE search_folders(id, depth, path) AS (
          SELECT f.id, 1, ARRAY[ #{
            folder_names.map do |name|
              db.literal(name)
            end.join(', ')
            } ]
            FROM folders f
            WHERE f.bucket_id = #{bucket.id} AND f.folder_name = #{db.literal(folder_names.first)} AND f.folder_id IS NULL
          UNION ALL
            SELECT f.id, sf.depth+1, sf.path
            FROM folders f, search_folders sf
            WHERE bucket_id = #{bucket.id}
              AND f.folder_name = sf.path[sf.depth+1]
              AND f.folder_id = sf.id
              AND sf.depth < #{folder_names.length + 1}
        )

        SELECT f.*
        FROM search_folders AS sf JOIN folders AS f
        ON f.id=sf.id ORDER BY sf.depth ASC;
      )

      # Find the set of folders that consecutively link as a chain from the root.
      folders = Metis::Folder.with_sql(query).all.inject([]) do |parents, folder|
        break parents unless folder.folder == parents.last
        parents << folder
      end

      if folders.length != folder_names.length
        return [] unless allow_partial_match
      end

      return folders
    end

    def self.exists?(folder_name, bucket, parent_folder)
      Metis::Folder.where(folder_name: folder_name, bucket: bucket, folder_id: parent_folder ? parent_folder.id : nil).count > 0
    end

    def parent_folders
      query =  %Q(
          WITH RECURSIVE parent_folders(id, folder_id, depth) AS (
            SELECT f.id, f.folder_id, 1
              FROM folders f
              WHERE f.bucket_id = #{bucket.id} AND f.id = #{db.literal(folder_id)}
            UNION ALL
              SELECT f.id, f.folder_id, pf.depth+1
              FROM folders f, parent_folders pf
              WHERE f.bucket_id = #{bucket.id} AND f.id = pf.folder_id AND depth < 256
          )
          SELECT f.*
          FROM parent_folders AS pf JOIN folders AS f
          ON f.id=pf.id ORDER BY pf.depth DESC;
        )
      Metis::Folder.with_sql(
        query
      ).all
    end

    def child_folders
      query =  %Q(
          WITH RECURSIVE child_folders(id, folder_id, depth) AS (
            SELECT f.id, f.folder_id, 1
              FROM folders f
              WHERE f.bucket_id = #{bucket.id} AND f.folder_id = #{db.literal(id)}
            UNION ALL
              SELECT f.id, f.folder_id, cf.depth+1
              FROM folders f, child_folders cf
              WHERE f.bucket_id = #{bucket.id} AND f.folder_id = cf.id AND depth < 256
          )
          SELECT f.*
          FROM child_folders AS cf JOIN folders AS f
          ON f.id=cf.id ORDER BY cf.depth ASC;
        )

      Metis::Folder.with_sql(
        query
      ).all
    end

    def read_only?
      read_only
    end

    def root_folder?
      folder_id == nil
    end

    def size
      Metis::File.total_size(
        folder: [ self ] + child_folders
      )
    end

    def folder_path
      parent_folders.map(&:folder_name) + [ folder_name ]
    end

    def can_remove?
      !read_only? && files.empty? && folders.empty?
    end

    def rename!(new_folder, new_folder_name, user=nil)
      new_params = {
        folder: new_folder,
        folder_name: new_folder_name
      }

      new_params[:author] = Metis::File.author(user) if user

      update(**new_params)
      refresh
    end

    def update_bucket_and_rename!(new_folder, new_folder_name, new_bucket, user=nil)
      new_params = {
        folder: new_folder,
        folder_name: new_folder_name,
        bucket: new_bucket,
      }

      new_params[:author] = Metis::File.author(user) if user

      update(**new_params)

      # Need to recursively update all sub-folders and files
      files.each { |file|
        file.update_bucket!(new_bucket)
      }
      folders.each { |folder|
        folder.update_bucket_and_rename!(self, folder.folder_name, new_bucket)
      }

      refresh
    end

    def copy(new_parent_folder: nil, new_bucket: nil, user:)
      raise CopyError.new("Cannot copy into child folder.") if new_parent_folder && child_folders.include?(new_parent_folder)
      raise CopyError.new("Cannot copy into itself.") if new_parent_folder && new_parent_folder == self

      new_folder = Metis::Folder.create(
        folder_name: folder_name,
        bucket_id: new_bucket.nil? ? bucket.id : new_bucket.id,
        folder_id: new_parent_folder&.id,
        project_name: project_name,
        author: Metis::File.author(user),
      )

      # Need to recursively update all sub-folders and files
      files.each { |file|
        copy_file(new_folder, file, user)
      }
      folders.each { |folder|
        folder.copy(new_parent_folder: new_folder, new_bucket: new_bucket, user: user)
      }
    end

    def remove_contents!
      # remove child files
      Metis::File.where(folder_id: [ id ] + child_folders.map(&:id)).delete

      # remove child folders
      Metis::Folder.where(id: child_folders.map(&:id)).delete
    end

    def remove!
      delete
    end

    def protect!
      update(read_only: true)
    end

    def unprotect!
      update(read_only: false)
    end

    def to_hash(with_path=true)
      my_hash = {
        folder_name: folder_name,
        bucket_name: bucket.name,
        folder_path: with_path ? ::File.join(folder_path) : nil,
        project_name: project_name,
        read_only: read_only,
        updated_at: updated_at.iso8601,
        created_at: created_at.iso8601,
        author: author,
        id: id
      }

      my_hash.delete(:folder_path) if !with_path

      return my_hash
    end

    private

    def copy_file(new_folder, file, user)
      revision = Metis::CopyRevision.new({
        source: Metis::Path.path_from_parts(
          project_name,
          bucket.name,
          ::File.join(folder_path, file.file_name)
        ),
        dest: Metis::Path.path_from_parts(
          project_name,
          new_folder.bucket.name,
          ::File.join(new_folder.folder_path, file.file_name)
        ),
        user: user
      })

      revision.set_bucket(revision.source, [bucket])
      revision.set_bucket(revision.dest, [new_folder.bucket])

      revision.set_folder(
        revision.source,
        [self])
      revision.set_folder(
        revision.dest,
        [new_folder])

      revision.validate

      raise Etna::BadRequest, revision.errors.join("\n") unless revision.valid?

      revision.revise!
    end
  end
end
