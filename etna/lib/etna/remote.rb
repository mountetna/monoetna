require "net/ssh"
require "net/scp"
require "tempfile"

module Etna
  class RemoteSSH
    class RemoteSSHError < Exception
    end

    def initialize(host:, username:, password: nil, port: 22, root:, **args)
      @username = username
      @password = password
      @host = host
      @port = port
      @root = root
    end

    def ssh
      @ssh ||= Net::SSH.start(@host, @username, password: @password)
    end

    def mkdir_p(dir)
      output = ssh.exec!("mkdir -p #{dir}")

      raise RemoteSSHError.new("Unable to mkdir -p, #{output}") unless 0 == output.exitstatus
    end

    def lftp_get(username:, password:, host:, remote_filename:, &block)
      full_local_path = ::File.join(@root, host, remote_filename)
      full_local_dir = ::File.dirname(full_local_path)
      mkdir_p(full_local_dir)

      cmd = "lftp sftp://#{username}:#{password}@#{host}  -e \"get #{remote_filename} -o #{full_local_path}; bye\""

      output = ssh.exec!(cmd)
      raise RemoteSSHError.new("LFTP get failure: #{output}") unless 0 == output.exitstatus
      yield remote_filename if block_given?
    end

    def file_upload(remote_path, content)
      begin
        Tempfile.create do |temp_file|
          temp_file.binmode
          temp_file.write(content)
          temp_file.flush
          ssh.scp.upload!(temp_file.path, remote_path)
        end
      rescue StandardError => e
        raise RemoteSSHError.new("File upload failed: #{e.message}")
      end
    end
  end
end
