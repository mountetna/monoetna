require 'digest'
require 'fileutils'
require 'securerandom'

class Vulcan
  class Storage
    attr_accessor :data_root, :data_host

    def initialize(
        data_host: Vulcan.instance.config(:vulcan)[:host],
        data_root: Vulcan.instance.config(:data_folder),
        advisory_lock_file: Vulcan.instance.config(:advisory_lock_file)
    )
      @data_host = data_host
      @data_root = data_root
      @advisory_lock_file = advisory_lock_file
    end

    # Raised whenever two processes attempt to build and write the same build output
    class ContentionError < StandardError
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

    def run_cell_for(buildable, uniq_run_id: SecureRandom.uuid)
      NonSharedTmpRunCell.new(
          project_name: buildable.project_name, build_cell_hash: buildable.cell_hash,
          output_names: buildable.build_outputs.keys, uniq_run_id: uniq_run_id
      )
    end

    # Provides a context for which to prepare a temporary, non shared run cell to execute some calculation into.
    # Cleans up the run cell by the end of the block.  Does NOT retain the storage lock, as the run_cell is
    # intended to be non shared.
    def with_run_cell_build(run_cell:)
      build_dir = cell_data_path(project_name: run_cell.project_name, cell_hash: run_cell.cell_hash, prefix: 'tmp')

      unless ::File.exists?(build_dir)
        ::FileUtils.mkdir_p(build_dir)
      end

      output_storage_files = run_cell.build_outputs.values.map do |sf|
        sf = sf.with_prefix('tmp')
        ::FileUtils.touch(sf.data_path(self))
        sf
      end

      yield build_dir, output_storage_files
    ensure
      ::FileUtils.rm_rf(build_dir) if ::File.exists?(build_dir)
    end

    # Attempts to install sthe given build_dir, usually obtained via with_run_cell_build, into the final
    # location for the given buildable.  A lock is obtained to ensure an atomic transaction, raising a ContentionError
    # when multiple concurrent writes contend for access.  Graceful handlers should catch this error and consider wether
    # there is useful action by the system (reporting?) but likely not bother an end user, as it is an expected result
    # of system behavior.
    def install_build_output(buildable:, build_dir:)
      outputs_dir = cell_data_path(project_name: buildable.project_name, cell_hash: buildable.cell_hash)

      with_lock do
        if ::File.exists?(outputs_dir)
          raise ContentionError.new("Contention writing #{outputs_dir}, another process has already built it.")
        end

        ::FileUtils.mkdir_p(::File.dirname(outputs_dir))
        ::FileUtils.mv(build_dir, outputs_dir)
      end
    end

    # Advisory only level lock that will only be respected by other magma processes, not any other.
    # Supports re-entry per thread via thread locals
    # In the future this should be replaced with a proper distributed synchronization primitive, maybe redis or
    # NFS 4's close to open consistency backing a 2 pass polling lock.
    def with_lock(&block)
      if Thread.current[:vulcan_storage_lock_acquired]
        yield
      else
        lock_file = ::File.open(@advisory_lock_file, File::CREAT | File::RDWR)
        lock_file.flock(::File::LOCK_EX)

        begin
          Thread.current[:vulcan_storage_lock_acquired] = true
          yield
        ensure
          Thread.current[:vulcan_storage_lock_acquired] = false
          lock_file.flock(::File::LOCK_UN)
        end
      end
    end

    # An input or output file contained within a cell.
    class StorageFile
      # The logical name of a file is either the same as the data_filename (in the case of output files), or an alias
      # name given as the input's name (in that case, the source of the input file is linked from its original data_filename
      # to the logical_name).
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

      def with_prefix(prefix)
        StorageFile.new(
            project_name: project_name, cell_hash: cell_hash,
            data_filename: data_filename, logical_name: logical_name,
            prefix: prefix)
      end

      def data_path(storage)
        storage.data_path(project_name: project_name, cell_hash: cell_hash, data_filename: data_filename, prefix: prefix)
      end

      def as_input(logical_name)
        StorageFile.new(
            project_name: project_name, cell_hash: cell_hash,
            data_filename: data_filename, logical_name: logical_name)
      end

      def to_archimedes_storage_file(storage)
        path = data_path(storage)
        if (host_dir_map = ENV['HOST_DIR_MAP'])
          container_path, host_dir = host_dir_map.split('=', 2)
          path.sub!(/^#{container_path}/, host_dir)
        end

        "#{logical_name}:#{path}"
      end
    end

    # While StorageFiles and BuildTargets are shared contexts that must be acted on atomically and carefully,
    # run cells are non shared, expected to be completely unique to a given scheduled request, and are completely
    # ephemeral in the strictest sense (could be deleted anytime).  They provide the context of a specific cell's
    # run which is, if the cell completes successful, installed into the final BuildTarget via install_build_output
    class NonSharedTmpRunCell
      attr_reader :project_name, :build_cell_hash, :output_names, :uniq_run_id

      def initialize(project_name:, build_cell_hash:, output_names:, uniq_run_id:)
        @project_name = project_name
        @build_cell_hash = build_cell_hash
        @output_names = output_names
        @uniq_run_id = uniq_run_id
      end

      def cell_hash
        @cell_hash ||= Storage.hash_json_obj({
            build_cell_hash: build_cell_hash,
            uniq_run_id: uniq_run_id
        })
      end

      def build_outputs
        @build_outputs ||= output_names.map do |name|
          [
              name,
              StorageFile.new(
                  project_name: project_name, cell_hash: cell_hash,
                  data_filename: name, logical_name: name, prefix: 'tmp'
              )
          ]
        end.to_h
      end
    end

    # A canonical data result of a cell.
    class BuildTarget
      attr_reader :script, :project_name, :input_files

      def initialize(project_name:, session_key:, input_files:, output_filenames:, script:)
        unless input_files.map(&:project_name).all? { |v| v == project_name }
          raise "input files are mixed across projects, they must all belong to #{project_name}"
        end

        input_files = input_files.dup.sort_by { |input_file| [input_file.cell_hash, input_file.data_filename, input_file.logical_name] }
        output_filenames = output_filenames.dup.sort

        @project_name = project_name
        @session_key = session_key
        @input_files = input_files
        @output_filenames = output_filenames
        @script = script
      end

      def cell_hash
        @cell_hash ||= Storage.hash_json_obj({
            project_name: @project_name,
            input_files: @input_files,
            output_filenames: @output_filenames,
            session_key: @session_key,
            script: @script,
        })
      end

      def build_outputs
        @build_outputs ||= @output_filenames.map do |output_filename|
          [
              output_filename,
              StorageFile.new(
                  project_name: @project_name, cell_hash: cell_hash,
                  data_filename: output_filename, logical_name: output_filename
              ),
          ]
        end.to_h
      end

      def take_as_input(output_name, input_name)
        build_outputs[output_name]&.as_input(input_name)
      end

      def is_buildable?(storage)
        input_files.all? { |input_file| ::File.exists?(input_file.data_path(storage)) }
      end

      def is_built?(storage)
        build_outputs.values.all? { |output_file| ::File.exists?(output_file.data_path(storage)) }
      end

      def should_build?(storage)
        is_buildable?(storage) && !is_built?(storage)
      end
    end

    class MaterialSource
      attr_reader :project_name

      def initialize(project_name:, session_key:, material_reference:)
        @project_name = project_name
        @session_key = session_key
        @digest = Storage.hash_json_obj(material_reference)
      end

      def cell_hash
        @cell_hash ||= Storage.hash_json_obj({
            project_name: @project_name,
            output_filenames: ['material.bin'],
            session_key: @session_key,
            digest: @digest,
        })
      end

      def build_outputs
        @build_outputs ||= {
            'material.bin' =>
                StorageFile.new(
                    project_name: @project_name,
                    cell_hash: cell_hash,
                    data_filename: 'material.bin',
                    logical_name: 'material.bin',
                )
        }
      end

      def build_output
        build_outputs.values.first
      end

      def take_as_input(input_name)
        build_output&.as_input(input_name)
      end
    end

    def self.hash_json_obj(obj)
      parts = []
      q = [['root', obj]]

      until q.empty?
        k, val = q.shift
        val = val.as_json if val.respond_to?(:as_json)

        if val.is_a?(Hash)
          val.keys.sort.each do |subk|
            q << ["#{k}.#{subk}", val[subk]] unless val[subk].nil?
          end
        elsif val.is_a?(Array)
          val.each_with_index { |val, i| q << ["#{k}[#{i}]", val] }
        elsif val.nil? || val.is_a?(Numeric) || val.is_a?(String) || val.is_a?(Symbol) || val.is_a?(TrueClass) || val.is_a?(FalseClass)
          parts << "#{k}=#{JSON.dump(val)}"
        else
          raise "#{val.class.name} is not json serializable within #{obj}"
        end
      end

      Digest::SHA1.hexdigest(parts.join('&'))
    end
  end
end