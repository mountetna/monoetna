require 'yaml'
require 'fileutils'
require 'open3'
require 'securerandom'
require 'concurrent-ruby'
require 'net/sftp'
require 'net/ssh'

module Etna
  # A class that encapsulates opening / reading file system entries that abstracts normal file access in order
  # to make stubbing, substituting, and testing easier.
  class Filesystem
    def with_writeable(dest, opts = 'w', size_hint: nil, &block)
      raise "with_writeable not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.open(dest, opts, &block)
    end

    def ls(dir)
      raise "ls not supported by #{self.class.name}" unless self.class == Filesystem
      ::Dir.entries(dir).select { |p| !p.start_with?('.') }.map do |path|
        ::File.file?(::File.join(dir, path)) ? [:file, path] : [:dir, path]
      end
    end

    def with_readable(src, opts = 'r', &block)
      raise "with_readable not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.open(src, opts, &block)
    end

    def mkdir_p(dir)
      raise "mkdir_p not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.mkdir_p(dir)
    end

    def rm_rf(dir)
      raise "rm_rf not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.rm_rf(dir)
    end

    def tmpdir
      raise "tmpdir not supported by #{self.class.name}" unless self.class == Filesystem
      ::Dir.mktmpdir
    end

    def exist?(src)
      raise "exist? not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.exist?(src)
    end

    def mv(src, dest)
      raise "mv not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.mv(src, dest)
    end

    def stat(src)
      raise "stat not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.stat(src)
    end

    class EmptyIO < StringIO
      def write(*args)
        # Do nothing -- always leave empty
      end
    end

    module WithPipeConsumer
      def mkio(file, opts, size_hint: nil, &block)
        rd, wd = IO.pipe

        pid = spawn(*mkcommand(rd, wd, file, opts, size_hint: size_hint))
        q = Queue.new

        closer = Thread.new do
          _, status = Process.wait2 pid
          q << status
        end

        begin
          if opts.include?('w')
            rd.close
            yield wd
            wd.close
          else
            wd.close
            yield rd
            rd.close
          end

          closer.join
        rescue => e
          wd.close
          rd.close
          Process.kill("HUP", pid)
          raise e
        end

        status = q.pop
        raise IOError.new("Failed to run external process, got status code #{status}") unless status.success?
      end
    end

    class AsperaCliFilesystem < Filesystem
      include WithPipeConsumer

      def initialize(ascli_bin:, ascp_bin:, host:, username:, password: nil, key_file: nil, port: 33001)
        @ascli_bin = ascli_bin
        @ascp_bin = ascp_bin
        @username = username
        @password = password
        @key_file = key_file
        @host = host
        @port = port

        @config_file = File.join(Dir.mktmpdir, "config.yml")
        config = {}
        config["config"] = {"version" => `#{ascli_bin} --version`.chomp}
        config["default"] = {"server" => "clifilesystem"}
        server_config = config["clifilesystem"] = {
            "url" => "ssh://#{host}:#{port}",
            "username" => username,
            "ssh_options" => {append_all_supported_algorithms: true},
        }

        if password
          server_config["password"] = password
        elsif key_file
          server_config["ssh_keys"] = key_file
        else
          raise "One of password or key_file must be provided"
        end

        ::File.open(@config_file, "w") do |file|
          file.write(config.to_yaml)
        end
      end

      def run_ascli_cmd(cmd, *opts)
        output, status = Open3.capture2(@ascli_bin, "server", cmd, *opts, "--format=json", "--config=#{@config_file}")

        if status.success?
          return JSON.parse(output)
        end

        nil
      end

      def with_writeable(dest, opts = 'w', size_hint: nil, &block)
        mkio(dest, opts, size_hint: size_hint, &block)
      end

      def with_readable(src, opts = 'r', &block)
        mkio(src, opts, &block)
      end

      def mkdir_p(dir)
        raise "Failed to mkdir #{dir}" unless run_ascli_cmd("mkdir", dir)
      end

      def rm_rf(dir)
        raise "Failed to rm_rf #{dir}" unless run_ascli_cmd("rm", dir)
      end

      def tmpdir
        tmpdir = "/Upload/Temp/#{SecureRandom.hex}"
        mkdir_p(tmpdir)
        tmpdir
      end

      def exist?(src)
        !run_ascli_cmd("ls", src).nil?
      end

      def mv(src, dest)
        raise "Failed to mv #{src} to #{dest}" unless run_ascli_cmd("mv", src, dest)
      end

      def mkcommand(rd, wd, file, opts, size_hint: nil)
        env = {}
        cmd = [env, @ascp_bin]

        if @password
          env['ASPERA_SCP_PASS'] = @password
        else
          cmd << "-i"
          cmd << @key_file
        end

        cmd << "-P"
        cmd << @port.to_s

        remote_path = file
        # https://download.asperasoft.com/download/docs/entsrv/3.9.1/es_admin_linux/webhelp/index.html#dita/stdio_2.html
        local_path = "stdio://"
        if size_hint
          local_path += "/?#{size_hint}"
        end

        if opts.include?('r')
          cmd << '--mode=recv'
          cmd << "--host=#{@host}"
          cmd << "--user=#{@username}"
          cmd << remote_path
          cmd << local_path

          cmd << {out: wd}
        elsif opts.include?('w')
          cmd << '--mode=send'
          cmd << "--host=#{@host}"
          cmd << "--user=#{@username}"
          cmd << local_path
          cmd << remote_path

          cmd << {in: rd}
        end

        cmd
      end
    end

    # Genentech's aspera deployment doesn't support modern commands, unfortunately...
    class GeneAsperaCliFilesystem < AsperaCliFilesystem
      def mkdir_p(dest)
        # Pass through -- this file system creates containing directories by default, womp womp.
      end

      def mkcommand(rd, wd, file, opts, size_hint: nil)
        if opts.include?('w')
          super.map do |e|
            if e.instance_of?(String) && e.start_with?("stdio://")
              "stdio-tar://"
            elsif e == file
              ::File.dirname(file)
            else
              e
            end
          end.insert(2, "-d")
        else
          super
        end
      end

      def with_writeable(dest, opts = 'w', size_hint: nil, &block)
        raise "#{self.class.name} requires size_hint in with_writeable" if size_hint.nil?

        super do |io|
          io.write("File: #{::File.basename(dest)}\n")
          io.write("Size: #{size_hint}\n")
          io.write("\n")
          yield io
        end
      end
    end

    class Metis < Filesystem
      def initialize(metis_client:, project_name:, bucket_name:, root: '/', uuid: SecureRandom.uuid)
        @metis_client = metis_client
        @project_name = project_name
        @bucket_name = bucket_name
        @root = root
        @metis_uid = uuid
      end

      def metis_path_of(path)
        joined = ::File.join(@root, path)
        joined[0] == "/" ? joined.slice(1..-1) : joined
      end

      def create_upload_workflow
        Etna::Clients::Metis::MetisUploadWorkflow.new(metis_client: @metis_client, metis_uid: @metis_uid, project_name: @project_name, bucket_name: @bucket_name, max_attempts: 3)
      end

      def with_hot_pipe(opts, receiver, *args, &block)
        rp, wp = IO.pipe
        begin
          executor = Concurrent::SingleThreadExecutor.new(fallback_policy: :abort)
          begin
            if opts.include?('w')
              future = Concurrent::Promises.future_on(executor) do
                self.send(receiver, rp, *args)
              rescue => e
                Etna::Application.instance.logger.log_error(e)
                raise e
              ensure
                rp.close
              end

              yield wp
            else
              future = Concurrent::Promises.future_on(executor) do
                self.send(receiver, wp, *args)
              rescue => e
                Etna::Application.instance.logger.log_error(e)
                raise e
              ensure
                wp.close
              end

              yield rp
            end

            future.wait!
          ensure
            executor.shutdown
            executor.kill unless executor.wait_for_termination(5)
          end
        ensure
          rp.close
          wp.close
        end
      end

      def do_streaming_upload(rp, dest, size_hint)
        streaming_upload = Etna::Clients::Metis::MetisUploadWorkflow::StreamingIOUpload.new(readable_io: rp, size_hint: size_hint)
        create_upload_workflow.do_upload(
            streaming_upload,
            metis_path_of(dest)
        )
      end

      def with_writeable(dest, opts = 'w', size_hint: nil, &block)
        self.with_hot_pipe(opts, :do_streaming_upload, dest, size_hint) do |wp|
          yield wp
        end
      end

      def create_download_workflow
        Etna::Clients::Metis::MetisDownloadWorkflow.new(metis_client: @metis_client, project_name: @project_name, bucket_name: @bucket_name, max_attempts: 3)
      end

      def do_streaming_download(wp, metis_file)
        create_download_workflow.do_download(wp, metis_file)
      end

      def with_readable(src, opts = 'r', &block)
        metis_file = list_metis_directory(::File.dirname(src)).files.all.find { |f| f.file_name == ::File.basename(src) }
        raise "Metis file at #{@project_name}/#{@bucket_name}/#{@root}/#{src} not found.  No such file" if metis_file.nil?

        self.with_hot_pipe(opts, :do_streaming_download, metis_file) do |rp|
          yield rp
        end
      end

      def list_metis_directory(path)
        @metis_client.list_folder(Etna::Clients::Metis::ListFolderRequest.new(project_name: @project_name, bucket_name: @bucket_name, folder_path: metis_path_of(path)))
      end

      def mkdir_p(dir)
        create_folder_request = Etna::Clients::Metis::CreateFolderRequest.new(
            project_name: @project_name,
            bucket_name: @bucket_name,
            folder_path: metis_path_of(dir),
        )
        @metis_client.create_folder(create_folder_request)
      end

      def ls(dir)
        response = list_metis_directory(::File.dirname(dir))
        response.files.map { |f| [:file, f.file_name] } + response.folders.map { |f| [:folder, f.folder_name] }
      end

      def exist?(src)
        begin
          response = list_metis_directory(::File.dirname(src))
        rescue Etna::Error => e
          if e.status == 404
            return false
          elsif e.message =~ /Invalid folder/
            return false
          end

          raise e
        end

        response.files.all.any? { |f| f.file_name == ::File.basename(src) } ||
            response.folders.all.any? { |f| f.folder_name == ::File.basename(src) }
      end
    end

    class SftpFilesystem < Filesystem
      def initialize(host:, username:, password: nil, port: 22, **args)
        @username = username
        @password = password
        @host = host
        @port = port
      end

      def ssh
        @ssh ||= Net::SSH.start(@host, @username, password: @password)
      end

      def sftp
        @sftp ||= begin
          conn = Net::SFTP::Session.new(ssh)
          conn.loop { conn.opening? }

          conn
        end
      end

      def with_readable(src, opts = 'r', &block)
        sftp.file.open(src, opts, &block)
      end

      def ls(dir)
        sftp.dir.entries(dir)
      end

      def exist?(src)
        begin
          sftp.file.open(src)
        rescue Net::SFTP::StatusException
          return false
        end
        return true
      end

      def stat(src)
        sftp.file.open(src).stat
      end
    end

    class RemoteSSH < Filesystem
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
    end

    class Mock < Filesystem
      class MockStat
        def initialize
        end

        def size
          0
        end
      end

      def initialize(&new_io)
        @files = {}
        @dirs = {}
        @new_io = new_io
      end

      def mkio(file, opts)
        if @new_io.nil?
          StringIO.new
        else
          @new_io.call(file, opts)
        end
      end

      def with_writeable(dest, opts = 'w', size_hint: nil, &block)
        if @dirs.include?(dest)
          raise IOError.new("Path #{dest} is a directory")
        end

        dir, file = File.split(dest)
        @dirs[dir] ||= Set.new
        @dirs[dir].add(file)
        yield (@files[dest] = mkio(dest, opts))
      end

      def mkdir_p(dest)
        while !@dirs.include?(dest)
          @dirs[dest] = Set.new
          break if dest == "." || dest == "/"
          dest, _ = File.split(dest)
        end
      end

      def mv(src, dest)
        if exist?(dest)
          raise "#{dest} already exists, cannot move"
        end

        if @dirs.include?(src)
          @dirs[dest] = @dirs.delete(src)
        elsif @files.include?(src)
          @files[dest] = @files.delete(src)
        else
          raise "#{src} does not exist, cannot move"
        end
      end

      def tmpdir
        require 'securerandom'
        "/tmp-#{SecureRandom::uuid}"
      end

      def with_readable(src, opts = 'r', &block)
        if @dirs.include?(src)
          raise IOError.new("Path #{src} is a directory")
        end

        if !@files.include?(src)
          raise IOError.new("Path #{src} does not exist")
        end

        @files[src].rewind
        yield @files[src]
      end

      def exist?(src)
        @files.include?(src) || @dirs.include?(src)
      end

      def stat(src)
        @files[src].respond_to?(:stat) ? @files[src].stat : MockStat.new
      end
    end
  end
end