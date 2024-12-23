require 'net/sftp'

class SFTPClient
  def initialize(host, username, password, port: 22)
    @host = host
    @username = username
    @password = password
    @port = port
    begin
      Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
        sftp.dir.entries('.')
      end
    rescue StandardError => e
      raise "Failed to establish SFTP connection: #{e.message}"
    end
  end

  # Uploads a local file to the remote SFTP server
  #
  # @param local_path [String] the path to the local file
  # @param remote_path [String] the destination path on the remote server
  def upload(local_path, remote_path)
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      sftp.upload!(local_path, remote_path)
      puts "Uploaded #{local_path} to #{remote_path}"
    end
  rescue StandardError => e
    puts "Upload failed: #{e.message}"
  end

  # Downloads a file from the remote SFTP server to the local machine
  #
  # @param remote_path [String] the path to the remote file
  # @param local_path [String] the destination path on the local machine
  def download(remote_path, local_path)
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      sftp.download!(remote_path, local_path)
      puts "Downloaded #{remote_path} to #{local_path}"
    end
  rescue StandardError => e
    puts "Download failed: #{e.message}"
  end

    # Recursively searches for files in remote directories that match a given regex pattern
    #
    # @param remote_dir [String] the remote directory path to start search from
    # @param pattern [Regexp] the regex pattern to match file names
    # @param start_time [Time] the start of the modified time range
    # @param end_time [Time] the end of the modified time range
    # @param ignore_dirs [Array<String>] list of directory names to ignore
    # @return [Array<String>] list of matching file paths
    def search_files(remote_dir, pattern, start_time, end_time, ignore_dirs=[])
      matching_files = []
      Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
        search_directory(sftp, remote_dir, pattern, start_time, end_time, ignore_dirs, matching_files)
      end
      matching_files
    rescue StandardError => e
      puts "Search failed: #{e.message}"
      []
    end

    private
    def search_directory(sftp, current_dir, pattern, start_time, end_time, ignore_dirs, matching_files)
      entries = sftp.dir.entries(current_dir)
      entries.each do |entry|
        next if entry.name == '.' || entry.name == '..'
        next if ignore_dirs.include?(entry.name)

        full_path = File.join(current_dir, entry.name)
        if entry.directory?
          search_directory(sftp, full_path, pattern, start_time, end_time, ignore_dirs, matching_files)
        elsif entry.name.match?(pattern) && entry.attributes.mtime >= start_time && entry.attributes.mtime <= end_time
          matching_files << { path: full_path, modified_time: entry.attributes.mtime }
        end
      end
    end

  # Lists files in a remote directory
  #
  # @param remote_dir [String] the remote directory path
  def list_files(remote_dir)
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      entries = sftp.dir.entries(remote_dir)
      entries.each { |entry| puts entry.name }
    end
  rescue StandardError => e
    puts "Listing files failed: #{e.message}"
  end
end