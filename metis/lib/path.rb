class Metis
    class Path

      attr_reader :path
      def initialize(path)
        @path = path
      end

      def self.filepath_match
        # Assumes Metis file paths are in the form
        #    metis://<project>/<bucket>/<folder path>/<file name>
        # Splitting the above produces
        #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]

        @filepath_match ||= %r!
        \A
          metis:\/\/
          (?<project_name>(\w+))
          /
          (?<bucket_name>(\w+))
          /
          (?<file_path>(#{Metis::File::FILENAME_MATCH.source}).*)
        \z
      !x
      end

      def self.root_folder_path_match
        # Assumes Metis root folder paths are in the form
        #    metis://<project>/<bucket>
        # Splitting the above produces
        #   ["metis", "", "<project>", "<bucket>"]

        @root_folder_path_match ||= %r!
        \A
          metis:\/\/
          (?<project_name>(\w+))
          /
          (?<bucket_name>(\w+))
          /?
        \z
      !x
      end

      def self.remove_double_slashes(path)
        return ::File.join(path.split("/")) unless path.start_with?("metis://")

        ::File.join("metis://", path.split("metis://").last.split("/"))
      end

      def self.path_from_parts(project_name, bucket_name, file_path)
        ::File.join("metis://", project_name, bucket_name, file_path.split("/"))
      end

      def project_name
        valid? ?
          path_regex[:project_name] :
          nil
      end

      def bucket_name
        valid? ? path_regex[:bucket_name] : nil
      end

      def bucket_matches?(bucket)
        project_name == bucket.project_name && bucket_name == bucket.name
      end

      def folder_matches?(folder)
        ::File.join(folder.folder_path) == folder_path &&
        folder.project_name == project_name
      end

      def folder_path
        folder_path, _ = Metis::File.path_parts(file_path)
        return folder_path
      end

      def file_path
        valid? && matches_filepath? ?
          Metis::Path.filepath_match.match(@path)[:file_path] :
          nil
      end

      def file_name
        _, file_name = Metis::File.path_parts(file_path)
        return file_name
      end

      def valid?
        matches_filepath? || matches_root_folder_path?
      end

      def matches_filepath?
        !!Metis::Path.filepath_match.match(@path)
      end

      def matches_root_folder_path?
        !!Metis::Path.root_folder_path_match.match(@path)
      end

      def path_regex
        matches_filepath? ?
          Metis::Path.filepath_match.match(@path) :
          Metis::Path.root_folder_path_match.match(@path)
      end
    end
  end
