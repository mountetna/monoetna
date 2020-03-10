class Metis
  class Assimilation
    def self.assimilate_folder(folder, file_path)
      self.new(
        file_path,
        folder.bucket,
        folder
      ).execute
    end

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
        Dir.glob(["#{@path}/*"]).each do |file_path|
          Metis::Assimilation.assimilate_folder(metis_folder, file_path)
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

      print_status(folder_exists? && 'folder exists, continuing')

      metis_folder = Metis::Folder.find_or_create(
        project_name: @bucket.project_name,
        folder_name: @name,
        bucket: @bucket,
        folder_id: @parent_folder&.id
      ) do |folder|
        folder.author = '|Metis'
      end

      return metis_folder
    end

    def ensure_file
      raise ArgumentError, 'Cannot write file over folder' if folder_exists?

      metis_file = Metis::File.find(
        project_name: @bucket.project_name,
        file_name: @name,
        folder_id: @parent_folder&.id,
        bucket: @bucket
      )

      if metis_file
        if ::File.size(@path) == metis_file.data_block.actual_size
          print_status('file exists, continuing')
        else
          print_status('file size mismatch, continuing')
        end
        return
      end

      print_status

      data_block = Metis::DataBlock.create_from(@name, @path, true)

      metis_file = Metis::File.create(
        project_name: @bucket.project_name,
        file_name: @name,
        folder_id: @parent_folder&.id,
        data_block: data_block,
        bucket: @bucket,
        author: '|Metis'
      )
    end

    def print_status(msg=nil)
      puts "#{@path} => /#{@parent_folder ? @parent_folder.folder_path.join('/') + '/' : nil }#{@name}#{msg ? " #{msg}" : nil}"
    end
  end
end
