class Metis
    class Path
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]

      attr_reader :path, :bucket
      def initialize(path)
        @path = path

        if valid?
          @bucket = Metis::Bucket.find(
            project_name: project_name,
            name: bucket_name
          )
        end
      end

      def filepath_match
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

      def self.path_from_parts(project_name, bucket_name, file_path)
        "metis://#{project_name}/#{bucket_name}/#{file_path}"
      end

      def project_name
        filepath_match.match(@path)[:project_name]
      end

      def bucket_name
        filepath_match.match(@path)[:bucket_name]
      end

      def folder
        Metis::Folder.from_path(
          @bucket, folder_path
        ).last
      end

      def folder_path
        folder_path, _ = Metis::File.path_parts(file_path)
        return folder_path
      end

      def file_path
        filepath_match.match(@path)[:file_path]
      end

      def file
        Metis::File.from_path(
          @bucket,
          file_path)
      end

      def file_name
        _, file_name = Metis::File.path_parts(file_path)
        return file_name
      end

      def valid?
        return false unless filepath_match.match(@path)
        return true
      end
    end
  end
