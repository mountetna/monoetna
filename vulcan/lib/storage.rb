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
    def with_build_transaction(buildable, &block)
      project_name = buildable.project_name
      ch = buildable.cell_hash
      build_dir = cell_data_path(project_name: project_name, cell_hash: ch, prefix: 'tmp')
      outputs_dir = cell_data_path(project_name: project_name, cell_hash: ch)

      if ::File.exists?(build_dir)
        ::FileUtils.rm_rf(build_dir)
      end

      ::FileUtils.mkdir_p(build_dir)

      output_storage_files = buildable.build_outputs.values.map do |sf|
        sf = sf.with_prefix('tmp')
        ::FileUtils.touch(sf.data_path(self))
        sf
      end

      result = yield output_storage_files

      ::FileUtils.mkdir_p(::File.dirname(outputs_dir))
      ::FileUtils.mv(build_dir, outputs_dir, force: true)

      result
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
          [output_filename,
              StorageFile.new(
                  project_name: @project_name, cell_hash: cell_hash,
                  data_filename: output_filename, logical_name: output_filename
              )]
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