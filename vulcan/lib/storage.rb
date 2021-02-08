require 'digest'
require 'fileutils'


class Vulcan
  class Storage
    attr_accessor :data_root, :data_host

    def initialize(
        data_host: Vulcan.instance.config(:vulcan)[:host],
        data_root: Vulcan.instance.config(:data_folder)
    )
      @data_host = data_host
      @data_root = data_root
    end

    # Data is namespaced per project, just like workflows and anything else, to make data access authorization much simpler.
    # You need project access to fetch any file under that path.
    def data_url(project_name:, cell_hash:, data_filename:)
      "#{data_host}/api/#{project_name}/data/#{cell_hash}/#{data_filename}"
    end

    def cell_data_path(project_name:, cell_hash:, prefix: "built")
      "#{data_root}/#{prefix}/#{project_name}/#{cell_hash}"
    end

    # Finalized location for built outputs
    def data_path(project_name:, cell_hash:, data_filename:, prefix: "built")
      "#{cell_data_path(project_name: project_name, cell_hash: cell_hash, prefix: prefix)}/#{data_filename}"
    end

    # Invoked to initialize temporary output files for a given cell definition, selecting a temporary data
    # location to place the tentative build files, and creating StorageFile instances to provide to the
    # build system.  After a successful invocation, places the temporary files into their true output locations.
    # When we move to a real scheduler, this 'finalization' logic will exist somewhere else (perhaps as a separate
    # endpoint for jobs to invoke themselves)
    def with_build_transaction(project_name:, ch:, output_filenames:, &block)
      build_dir = cell_data_path(project_name: project_name, cell_hash: ch, prefix: 'tmp')
      outputs_dir = cell_data_path(project_name: project_name, cell_hash: ch)

      if ::File.exists?(build_dir)
        ::FileUtils.rm_rf(build_dir)
      end

      ::FileUtils.mkdir_p(build_dir)

      output_storage_files = output_filenames.map do |filename|
        StorageFile.new(project_name: project_name, cell_hash: ch, data_filename: filename, logical_name: filename, prefix: 'tmp').tap do |sf|
          ::FileUtils.touch(sf.data_path(self))
        end
      end

      yield output_storage_files

      ::FileUtils.mkdir_p(::File.dirname(outputs_dir))
      ::FileUtils.mv(build_dir, outputs_dir, force: true)
    end

    class StorageFile
      # Logical name is simply a mapping between the canonical filename and the filename used inside of the cell context
      # For output files, these are the same
      # For input files, these are different
      attr_accessor :project_name, :cell_hash, :data_filename, :logical_name, :prefix

      def initialize(project_name:, cell_hash:, data_filename:, logical_name:, prefix: "built")
        @project_name = project_name
        @cell_hash = cell_hash
        @data_filename = data_filename
        @logical_name = logical_name
        @prefix = prefix
      end

      def as_json
        {
            project_name: @project_name,
            cell_hash: @cell_hash,
            data_filename: @data_filename,
            logical_name: @logical_name,
        }
      end

      def data_path(storage)
        storage.data_path(project_name: project_name, cell_hash: cell_hash, data_filename: data_filename, prefix: prefix)
      end

      def to_archimedes_storage_file(storage)
        "#{logical_name}:#{data_path(storage)}"
      end
    end

    def self.cell_hash(project_name:, input_files:, output_filenames:, session_key:, script: nil, raw_hash: nil)
      unless input_files.map(&:project_name).all? { |v| v == project_name }
        raise "input files are mixed across projects, they must all belong to #{project_name}"
      end

      input_files = input_files.dup.sort_by(&:logical_name)
      output_filenames = output_filenames.dup.sort

      hash_json_obj({
          project_name: project_name,
          input_files: input_files,
          output_filenames: output_filenames,
          session_key: session_key,
          script: script,
          raw_hash: raw_hash,
      })
    end

    def self.hash_json_obj(obj)
      parts = []
      q = [['root', obj]]

      until q.empty?
        k, val = q.shift
        val = val.as_json if val.respond_to?(:as_json)

        if val.is_a?(Hash)
          val.keys.sort.each do |subk|
            q << ["#{k}.#{subk}", val[subk]]
          end
        elsif val.is_a?(Array)
          val.each_with_index { |val, i| q << ["#{k}[#{i}]", val] }
        elsif val.nil? || val.is_a?(Numeric) || val.is_a?(String) || val.is_a?(TrueClass) || val.is_a?(FalseClass)
          parts << "#{k}=#{JSON.dump(val)}"
        else
          raise "#{val.class.name} is not json serializable"
        end
      end

      Digest::SHA1.hexdigest(parts.join('&'))
    end
  end
end