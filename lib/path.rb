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

      attr_reader :path
      def initialize(path)
        @path = path
      end

      def self.extract_bucket_name_from_path(path)
        return FILEPATH_MATCH.match(path)[:bucket_name]
      end

      def self.extract_file_path_from_path(path)
        return FILEPATH_MATCH.match(path)[:file_path]
      end

      def self.path_from_parts(project_name, bucket_name, file_path)
        "metis://#{project_name}/#{bucket_name}/#{file_path}"
      end

      def bucket_name
        Metis::Path.extract_bucket_name_from_path(@path)
      end

      def file_path
        Metis::Path.extract_file_path_from_path(@path)
      end

      def file
        Metis::File.from_path(
          get_bucket,
          file_path)
      end

      def valid?
        return false unless FILEPATH_MATCH.match(@path)
        return true
      end

      private

      def get_bucket
        Metis::Bucket.where(
          project_name: FILEPATH_MATCH.match(@path)[:project_name],
          name: FILEPATH_MATCH.match(@path)[:bucket_name]
        ).first
      end
    end
  end
