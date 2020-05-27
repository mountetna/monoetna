class Metis
    class Path
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]

      FILENAME_MATCH=/[^<>:;,?"*\|\/\x00-\x1f]+/x

      FILEPATH_MATCH=%r!
        \A
          metis:\/\/
          (?<project_name>(#{FILENAME_MATCH.source}))
          /
          (?<bucket_name>(#{FILENAME_MATCH.source}))
          /
          (?<file_path>(#{FILENAME_MATCH.source}).*)
        \z
      !x

      attr_reader :path, :bucket
      def initialize(path)
        @path = path

        if valid?
          @bucket = Metis::Bucket.find(
            project_name: FILEPATH_MATCH.match(@path)[:project_name],
            name: bucket_name
          )
        end
      end

      def self.path_from_parts(project_name, bucket_name, file_path)
        "metis://#{project_name}/#{bucket_name}/#{file_path}"
      end

      def project_name
        FILEPATH_MATCH.match(@path)[:project_name]
      end

      def bucket_name
        FILEPATH_MATCH.match(@path)[:bucket_name]
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
        FILEPATH_MATCH.match(@path)[:file_path]
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
        return false unless FILEPATH_MATCH.match(@path)
        return true
      end
    end
  end
