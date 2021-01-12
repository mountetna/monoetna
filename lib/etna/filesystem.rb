require 'yaml'
require 'fileutils'
require 'open3'

module Etna
  # A class that encapsulates opening / reading file system entries that abstracts normal file access in order
  # to make stubbing, substituting, and testing easier.
  class Filesystem
    def with_writeable(dest, opts = 'w', size_hint: nil, &block)
      ::File.open(dest, opts, &block)
    end

    def with_readable(src, opts = 'r', &block)
      ::File.open(src, opts, &block)
    end

    def mkdir_p(dir)
      ::FileUtils.mkdir_p(dir)
    end

    def rm_rf(dir)
      ::FileUtils.rm_rf(dir)
    end

    def tmpdir
      ::Dir.mktmpdir
    end

    def exist?(src)
      ::File.exist?(src)
    end

    def mv(src, dest)
      ::FileUtils.mv(src, dest)
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

          cmd << { out: wd }
        elsif opts.include?('w')
          cmd << '--mode=send'
          cmd << "--host=#{@host}"
          cmd << "--user=#{@username}"
          cmd << local_path
          cmd << remote_path

          cmd << { in: rd }
        end

        cmd
      end
    end

    class Mock < Filesystem
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
    end
  end
end