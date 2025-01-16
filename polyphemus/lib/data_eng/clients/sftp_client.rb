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

  def upload(local_path, remote_path)
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      sftp.upload!(local_path, remote_path)
      puts "Uploaded #{local_path} to #{remote_path}"
    end
  rescue StandardError => e
    puts "Upload failed: #{e.message}"
  end

  def download_as_stream(remote_path)
    io = StringIO.new
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      sftp.download!(remote_path, io)
    end
    io.rewind
    io
  rescue StandardError => e
    puts "Download as stream failed for #{remote_path}: #{e.message}"
    nil
  end

  # Recursively searches for files in remote directories that match a given regex pattern
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