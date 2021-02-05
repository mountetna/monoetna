require 'digest'
require 'fileutils'


class Vulcan
  module Storage
    # Data is namespaced per project, just like workflows and anything else, to make data access authorization much simpler.
    # You need project access to fetch any file under that path.
    def self.data_url(project_name:, cell_hash:, data_filename:)
      vulcan_host = Vulcan.instance.config(:vulcan)[:host]
      "#{vulcan_host}/api/#{project_name}/data/#{cell_hash}/#{data_filename}"
    end

    def self.cell_data_path(project_name:, cell_hash:, prefix: "built")
      data_path = Vulcan.instance.config(:data_path)
      "#{data_path}/#{prefix}/#{project_name}/#{cell_hash}"
    end

    # Finalized location for built outputs
    def self.data_path(project_name:, cell_hash:, data_filename:, prefix: "built")
      "#{cell_data_path(project_name: project_name, cell_hash: cell_hash, prefix: prefix)}/#{data_filename}"
    end

    # Invoked to initialize temporary output files for a given cell definition, selecting a temporary data
    # location to place the tentative build files, and creating StorageFile instances to provide to the
    # build system.  After a successful invocation, places the temporary files into their true output locations.
    # When we move to a real scheduler, this 'finalization' logic will exist somewhere else (perhaps as a separate
    # endpoint for jobs to invoke themselves)
    def self.with_output_transaction(project_name:, input_files:, output_filenames:, session_key:, script: nil, raw: nil, &block)
      ch = cell_hash(
          project_name: project_name, input_files: input_files,
          output_filenames: output_filenames, session_key: session_key,
          script: script, raw: raw
      )

      build_dir = cell_data_path(project_name: project_name, cell_hash: ch, prefix: 'tmp')
      outputs_dir = cell_data_path(project_name: project_name, cell_hash: ch)

      if ::File.exists?(build_dir)
        ::FileUtils.rm_rf(build_dir)
      end

      ::FileUtils.mkdir_p(build_dir)

      output_storage_files = output_filenames.map do |filename|
        StorageFile.new(project_name: project_name, cell_hash: ch, data_filename: filename, logical_name: filename).tap do |sf|
          ::File.touch(sf.data_path(prefix: 'tmp'))
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
      attr_accessor :project_name, :cell_hash, :data_filename, :logical_name

      def initialize(project_name:, cell_hash:, data_filename:, logical_name:)
        @project_name = project_name
        @cell_hash = cell_hash
        @data_filename = data_filename
        @logical_name = logical_name
      end

      def as_logical_key
        "#{project_name}/#{cell_hash}/#{data_filename}/#{logical_name}"
      end

      def data_path(**args)
        Storage.data_path(project_name: project_name, cell_hash: cell_hash, data_filename: data_filename, **args)
      end

      def to_archimedes_storage_file
        "#{logical_name}:#{data_path}"
      end
    end

    def self.cell_hash(project_name:, input_files:, output_filenames:, session_key:, script_or_raw_hash:)
      unless input_files.map(&:project_name).all? { |v| v == project_name }
        raise "input files are mixed across projects, they must all belong to #{project_name}"
      end

      keys = input_files.map(&:as_logical_key)
      keys.push(*output_filenames)
      keys << project_name
      keys << script_or_raw_hash
      keys << session_key
      # Sort so that input ordering does not change
      keys.sort!

      Digest::SHA1.hexdigest(keys.map { |k| URI.encode_www_form_component(k) }.join("&"))
    end
  end
end