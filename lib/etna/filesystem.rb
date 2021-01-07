module Etna
  # A class that encapsulates opening / reading file system entries that abstracts normal file access in order
  # to make stubbing, substituting, and testing easier.
  class Filesystem
    def with_writeable(dest, opts = 'w', &block)
      ::File.open(dest, opts, &block)
    end

    def with_readable(src, opts = 'r', &block)
      ::File.open(src, opts, &block)
    end

    def mkdir_p(dir)
      require 'fileutils'
      ::FileUtils.mkdir_p(dir)
    end

    def rm_rf(dir)
      require 'fileutils'
      FileUtils.rm_rf(dir)
    end

    def tmpdir
      ::Dir.tmpdir
    end

    def exist?(src)
      ::File.exist?(src)
    end

    class EmptyIO < StringIO
      def write(*args)
        # Do nothing -- always leave empty
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

      def with_writeable(dest, opts = 'w', &block)
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

      def tmpdir
        require 'securerandom'
        "/tmp-#{SecureRandom::uuid}"
      end

      def with_readable(src, opts = 'r', &block)
        if @dirs.include?(dest)
          raise IOError.new("Path #{dest} is a directory")
        end

        if !@files.include?(dest)
          raise IOError.new("Path #{dest} does not exist")
        end

        yield (@files[dest] ||= mkio(src, opts))
      end

      def exist?(src)
        @files.include?(src) || @dirs.include?(src)
      end
    end
  end
end