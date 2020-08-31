require_relative '../directed_graph'

module Etna
  module Codegen
    # Simple container class that collects inherited context, making deeply nested objects
    # easier to manage, while providing simple memoization of sub resources by key.
    class Context
      attr_reader :key

      def initialize(key)
        @key = {}
        @values = {}
        @key = prepare_key(key)
      end

      protected

      def sub_context(constructor, sub_key = nil, &block)
        if sub_key.nil?
          sub_key = constructor
          constructor = nil
        end

        key = prepare_key(sub_key)
        @values[sub_key] ||= begin
          if constructor.nil?
            result = block.call(key)
          else
            result = constructor.new(key)
          end

          notify(result) if self.respond_to? :notify
          result
        end
      end

      def prepare_key(key)
        raise "key must be a Hash, found #{key.class}" unless key.is_a? Hash
        overlap = (key.keys & @key.keys)
        raise "sub_resource keys overlap with parent: #{self} - #{overlap}" unless overlap.empty?
        @key.dup.update(key)
      end

      def lookup(sub_key)
        @values[sub_key]
      end
    end

    class Project < Context
      def initialize(package_name)
        super(project: self, package_name: @package_name)
        @build_path = ::Tempfile.mkdir
        @package_name = package_name
      end

      protected

      def build_file(file_path, &block)
        if file_path.start_with?('.') || file_path.start_with?('/')
          raise "file_path #{file_path} is not valid, cannot start with . or /"
        end

        sub_context(project_path: ProjectPath.new(self, file_path)) do |key|
          block.call(key)
        end
      end

      def files
        @values.filter { |k, v| k.include? :file }.map { |k, v| v }
      end
    end

    class ProjectPath
      attr_reader :project, :path

      def initialize(project, path)
        @project = project
        @path = path
      end

      def relative_path_to(other)
        adjusted_root = "."
        if @project != other.project
          if @project.package_name != other.project.package_name
            return [other.project.package_name] + other.path.split('/')
          end

          adjusted_root = Pathname.new(other.project.target_path).relative_path_from(@project.target_path).to_s
        end

        File.join(adjusted_root, Pathname.new(other.path).relative_path_from(@path).to_s).split('/')
      end
    end
  end
end
