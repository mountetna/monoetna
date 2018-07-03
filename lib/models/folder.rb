class Metis
  class Folder < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :bucket

    many_to_one :folder

    one_to_many :files

    plugin :rcte_tree, key: :folder_id

    def self.from_path(bucket, folder_path)
      return [] unless folder_path && bucket

      folder_names = folder_path.split(%r!/!)
      parent_folder_id = nil

      db = Metis.instance.db

      query = %Q(
        WITH RECURSIVE search_folders(id, depth, path) AS (
          SELECT id, 1, ARRAY[ #{
            folder_names.map do |name|
              db.literal(name)
            end.join(', ')
            } ]
            FROM folders
            WHERE bucket_id = #{bucket.id} AND folder_name = #{db.literal(folder_names.first)}
          UNION ALL
            SELECT f.id, sf.depth+1, sf.path
            FROM folders f, search_folders sf
            WHERE bucket_id = #{bucket.id} AND folder_name = sf.path[sf.depth+1] AND sf.depth < #{folder_names.length + 1}
        )
        SELECT f.*
        FROM search_folders AS sf JOIN folders AS f
        ON f.id=sf.id;
      )

      folders = Metis::Folder.with_sql(query).all

      root_folder = folders.first

      return [] unless !root_folder || root_folder.folder_id == nil

      return folders
    end

    def location
      ::File.expand_path(::File.join(
        bucket.location,
        folder_path.map do |folder_name|
          Metis::File.safe_file_name(folder_name)
        end
      ))
    end

    def read_only?
      read_only
    end

    def folder_path
      ancestors.map(&:folder_name) + [ folder_name ]
    end

    def has_directory?
      ::File.exists?(location) && ::File.directory?(location)
    end

    def can_remove?
      has_directory? && Dir.entries(location).size <= 2 && !read_only?
    end

    def remove!
      ::Dir.delete(location)
      delete
    end

    def protect!
      self.read_only = true
      self.save
    end

    def unprotect!
      self.read_only = false
      self.save
    end

    def create_actual_folder!
      FileUtils.mkdir_p(location)
    end

    def to_hash
      {
        folder_name: folder_name,
        bucket_name: bucket.name,
        project_name: project_name,
        updated_at: updated_at,
        created_at: created_at,
        author: author
      }
    end
  end
end
