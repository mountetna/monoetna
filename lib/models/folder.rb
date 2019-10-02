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

    class Assimilation
      def initialize(path, bucket, parent_folder=nil)
        @path = ::File.expand_path(path)
        @name = ::File.basename(@path)
        @bucket = bucket
        @parent_folder = parent_folder

        raise ArgumentError, "No such file #{@path}" unless ::File.exists?(@path)
        raise ArgumentError, "Invalid name #{@name}" unless Metis::File.valid_file_name?(@name)
      end

      def execute
        if folder?
          metis_folder = ensure_folder
          Dir.glob(["#{@path}/*"]).each do |file|
            metis_folder.assimilate(file)
          end
        else
          ensure_file
        end
      end

      private

      def folder?
        ::File.directory?(@path)
      end

      def folder_exists?
        Metis::Folder.exists?(@name, @bucket, @parent_folder)
      end

      def file_exists?
        Metis::File.exists?(@name, @bucket, @parent_folder)
      end

      def ensure_folder
        # create the folder
        raise ArgumentError, 'Cannot write folder over file' if file_exists?

        print_status(folder_exists? && 'exists, continuing')

        metis_folder = Metis::Folder.find_or_create(
          project_name: @bucket.project_name,
          folder_name: @name,
          bucket: @bucket,
          folder_id: @parent_folder&.id
        ) do |folder|
          folder.author = '|Metis'
        end

        metis_folder.create_actual_folder!

        return metis_folder
      end

      def ensure_file
        print_status

        raise ArgumentError, 'Cannot write file over folder' if folder_exists?

        metis_file = Metis::File.find(
          project_name: @bucket.project_name,
          file_name: @name,
          folder_id: @parent_folder&.id,
          bucket: @bucket
        )

        if metis_file
          if ::File.size(@path) == metis_file.actual_size
            print_status('file exists, continuing')
          else
            print_status('file size mismatch, continuing')
          end
          return
        end

        metis_file = Metis::File.create(
          project_name: @bucket.project_name,
          file_name: @name,
          folder_id: @parent_folder&.id,
          bucket: @bucket,
          author: '|Metis'
        )
        metis_file.set_file_data(@path, true)
      end

      def print_status(msg=nil)
        puts "#{@path} => /#{@parent_folder ? @parent_folder.folder_path.join('/') + '/' : nil }/#{@name}#{msg ? " #{msg}" : nil}"
      end
    end

    def assimilate(file_path)
      Metis::Folder::Assimilation.new(file_path, self.bucket, self).execute
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
