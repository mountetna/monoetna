class Metis
  class Folder < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :bucket

    many_to_one :folder

    one_to_many :folders

    one_to_many :files

    def self.from_path(bucket, folder_path)
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
            WHERE f.bucket_id = #{bucket.id} AND f.folder_name = #{db.literal(folder_names.first)}
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

      folders = Metis::Folder.with_sql(query).all

      # make sure we have found a complete path
      return [] unless !folders.empty? && folders.first.root_folder? && folders.length == folder_names.length

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

    def root_folder?
      folder_id == nil
    end

    def folder_path
      parent_folders.map(&:folder_name) + [ folder_name ]
    end

    def has_directory?
      ::File.exists?(location) && ::File.directory?(location)
    end

    def can_remove?
      has_directory? && Dir.entries(location).size <= 2 && !read_only?
    end

    def rename!(new_folder, new_folder_name)
      old_location = location

      update(folder: new_folder, folder_name: new_folder_name)
      refresh

      new_location = location
      FileUtils.mv(old_location, new_location)
    end

    def remove!
      ::Dir.delete(location)
      delete
    end

    def protect!
      update(read_only: true)
    end

    def unprotect!
      update(read_only: false)
    end

    def create_actual_folder!
      FileUtils.mkdir_p(location)
    end

    def self.assimilate(file_path, bucket, parent_folder=nil)
      folder_path = file_path = ::File.expand_path(file_path)
      folder_name = file_name = ::File.basename(file_path)

      raise ArgumentError unless ::File.exists?(file_path)
      raise ArgumentError if Metis::File.exists?(file_name, bucket, parent_folder) || Metis::Folder.exists?(folder_name, bucket, parent_folder)

      if ::File.directory?(folder_path)
        # create the folder
        metis_folder = Metis::Folder.create(
          project_name: bucket.project_name,
          folder_name: folder_name,
          bucket: bucket,
          folder: parent_folder,
          author: ''
        )
        metis_folder.create_actual_folder!

        puts ">#{folder_path} => #{parent_folder ? parent_folder.folder_path.join('/') : nil}#{folder_name}"

        Dir.glob(["#{folder_path}/*"]).each do |file|
          metis_folder.assimilate(file)
        end
      else
        puts ".#{file_path} => #{parent_folder ? parent_folder.folder_path.join('/') : nil}#{file_name}"
        metis_file = Metis::File.find_or_create(
          project_name: bucket.project_name,
          file_name: file_name,
          folder_id: parent_folder ? parent_folder.id : nil,
          bucket: bucket,
          author: ''
        )
        metis_file.set_file_data(file_path, true)
      end
    end

    def assimilate(file_path)
      Metis::Folder.assimilate(file_path, self.bucket, self)
    end

    def to_hash
      {
        folder_name: folder_name,
        bucket_name: bucket.name,
        folder_path: ::File.join(folder_path),
        project_name: project_name,
        read_only: read_only,
        updated_at: updated_at.iso8601,
        created_at: created_at.iso8601,
        author: author
      }
    end
  end
end
