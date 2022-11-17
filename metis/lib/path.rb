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

      def self.remove_double_slashes(path)
        ::File.join(path.split("/"))
      end

      def self.path_from_parts(project_name, bucket_name, file_path)
        ::File.join("metis://", project_name, bucket_name, file_path)
      end

      def project_name
        valid? ? Metis::Path.filepath_match.match(@path)[:project_name] : nil
      end

      def bucket_name
        valid? ? Metis::Path.filepath_match.match(@path)[:bucket_name] : nil
      end

      def bucket_matches?(bucket)
        project_name == bucket.project_name && bucket_name == bucket.name
      end

      def folder_matches?(folder)
        folder.folder_path.join('/') == folder_path &&
        folder.project_name == project_name
      end

      def folder_path
        folder_path, _ = Metis::File.path_parts(file_path)
        return folder_path
      end

      def file_path
        valid? ? Metis::Path.filepath_match.match(@path)[:file_path] : nil
      end

      def file_name
        _, file_name = Metis::File.path_parts(file_path)
        return file_name
      end

      def valid?
        !!Metis::Path.filepath_match.match(@path)
      end
    end
  end
