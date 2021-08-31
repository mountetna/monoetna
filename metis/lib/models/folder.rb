class Metis
  class Folder < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :bucket

    many_to_one :folder

    one_to_many :folders

    one_to_many :files

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

    def folder_path
      parent_folders.map(&:folder_name) + [ folder_name ]
    end

    def can_remove?
      !read_only? && files.empty?
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
  end
end
