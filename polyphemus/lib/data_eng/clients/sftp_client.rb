require 'net/sftp'

class SFTPClient
  def initialize(host, username, password, port: 22)
    @host = host
    @username = username
    @password = password
    @port = port
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

    # Searches for files in a remote directory that match a given regex pattern
  #
  # @param remote_dir [String] the remote directory path
  # @param pattern [Regexp] the regex pattern to match file names
  # @return [Array<String>] list of matching file names
  def search_files(remote_dir, pattern)
    matching_files = []
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      entries = sftp.dir.entries(remote_dir)
      entries.each do |entry|
        # Skip directories
        next if entry.directory?

        if entry.name.match?(pattern)
          matching_files << entry.name
          puts "Matched: #{entry.name}"
        end
      end
    end
    matching_files
  rescue StandardError => e
    puts "Search failed: #{e.message}"
    []
  end

  # Lists files in a remote directory
  #
  # @param remote_dir [String] the remote directory path
  def list_files(remote_dir)
    Net::SFTP.start(@host, @username, password: @password, port: @port) do |sftp|
      require 'pry'; binding.pry
      entries = sftp.dir.entries(remote_dir)
      entries.each { |entry| puts entry.name }
    end
  rescue StandardError => e
    puts "Listing files failed: #{e.message}"
  end
end