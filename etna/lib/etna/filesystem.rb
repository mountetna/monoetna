require 'yaml'
require 'fileutils'
require 'open3'
require 'securerandom'
require 'concurrent-ruby'
require 'curb'

module Etna
  # A class that encapsulates opening / reading file system entries that abstracts normal file access in order
  # to make stubbing, substituting, and testing easier.
  class Filesystem
    def with_writeable(dest, opts = 'w', size_hint: nil, &block)
      raise "with_writeable not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.open(dest, opts, &block)
    end

    class Error < StandardError
    end

    def ls(dir)
      raise Etna::Filesystem::Error, "ls not supported by #{self.class.name}" unless self.class == Filesystem
      ::Dir.entries(dir).select { |p| !p.start_with?('.') }.map do |path|
        ::File.file?(::File.join(dir, path)) ? [:file, path] : [:dir, path]
      end
    end

    def with_readable(src, opts = 'r', &block)
      raise Etna::Filesystem::Error, "with_readable not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.open(src, opts, &block)
    end

    def mkdir_p(dir)
      raise Etna::Filesystem::Error, "mkdir_p not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.mkdir_p(dir)
    end

    def rm_rf(dir)
      raise Etna::Filesystem::Error, "rm_rf not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.rm_rf(dir)
    end

    def tmpdir
      raise Etna::Filesystem::Error, "tmpdir not supported by #{self.class.name}" unless self.class == Filesystem
      ::Dir.mktmpdir
    end

    def exist?(src)
      raise Etna::Filesystem::Error, "exist? not supported by #{self.class.name}" unless self.class == Filesystem
      ::File.exist?(src)
    end

    def mv(src, dest)
      raise Etna::Filesystem::Error, "mv not supported by #{self.class.name}" unless self.class == Filesystem
      ::FileUtils.mv(src, dest)
    end

    def stat(src)
      raise Etna::Filesystem::Error, "stat not supported by #{self.class.name}" unless self.class == Filesystem
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

        cmd = mkcommand(rd, wd, file, opts, size_hint: size_hint)
        pid = spawn(*cmd)
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

    class Metis < Filesystem
      attr_reader :project_name, :bucket_name

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
      include WithPipeConsumer

      class SftpFile
        attr_reader :size, :name

        def initialize(metadata)
          @metadata_parts = metadata.split(" ")
          @size = @metadata_parts[4].to_i
          @perms = @metadata_parts.first
          @name = @metadata_parts[8]
        end
      end

      def initialize(host:, username:, password: nil, port: 22, **args)
        @username = username
        @password = password
        @host = host
        @port = port

        @dir_listings = {}
      end

      def url(src)
        "sftp://#{@host}/#{src}"
      end

      def authn
        "#{@username}:#{@password}"
      end

      def curl_cmd(path, opts=[])
        connection = Curl::Easy.new(url(path))
        connection.http_auth_types = :basic
        connection.username = @username
        connection.password = @password

        connection
      end

      def sftp_file_from_path(src)
        files = ls(::File.dirname(src)).split("\n").map do |listing|
          SftpFile.new(listing)
        end

        file = files.select do |file|
          file.name == ::File.basename(src)
        end

        raise "#{src} not found" if file.empty?

        file.first
      end

      def mkcommand(rd, wd, file, opts, size_hint: nil)
        env = {}
        cmd = [env, "curl"]

        cmd << "-u"
        cmd << authn
        cmd << "-o"
        cmd << "-"
        cmd << "-s"
        cmd << "-N"
        cmd << url(file)

        if opts.include?('r')
          cmd << {out: wd}
        end

        cmd
      end

      def with_readable(src, opts = 'r', &block)
        raise "#{src} does not exist" unless exist?(src)

        sftp_file = sftp_file_from_path(src)

        mkio(src, opts, size_hint: sftp_file.size, &block)
      end

      def exist?(src)
        files = ls(::File.dirname(src))
        files.include?(::File.basename(src))
      end

      def ls(dir)
        dir = dir + "/" unless "/" == dir[-1]  # Trailing / makes curl list directory

        return @dir_listings[dir] if @dir_listings.has_key?(dir)

        listing = ''
        connection = curl_cmd(dir)
        connection.on_body { |data| listing << data; data.size }
        connection.perform

        @dir_listings[dir] = listing

        listing
      end

      def stat(src)
        sftp_file_from_path(src)
      end
    end

    class Mock < Filesystem
      class MockStat
        def initialize(io)
          @io = io
        end

        def size
          @io.respond_to?(:length) ? @io.length : 0
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
        @files[src].respond_to?(:stat) ? @files[src].stat : MockStat.new(@files[src])
      end
    end
  end
end
