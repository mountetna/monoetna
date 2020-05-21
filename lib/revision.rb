require 'pry'
class Metis
  class Revision
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

    attr_reader :source, :dest
    def initialize(params)
      raise Etna::BadRequest, 'All revisions require "source" parameter' unless params.has_key? :source
      raise Etna::BadRequest, 'All revisions require "dest" parameter' unless params.has_key? :dest

      raise Etna::BadRequest, 'Invalid source path' unless valid_file_path?(params[:source])

      @source = params[:source]
      @dest = params[:dest]
      # TODO: Will have to handle DELETE differently,
      #       when params[:dest] will be nil
    end

    def self.create_from_parts(params)
      Metis::Revision.new({
        source: path_from_parts(
          params[:source][:project_name],
          params[:source][:bucket_name],
          params[:source][:file_path]),
        dest: path_from_parts(
          params[:dest][:project_name],
          params[:dest][:bucket_name],
          params[:dest][:file_path])
      })
    end

    def self.extract_bucket_name_from_path(path)
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]
      # Should this be in a central gem, like etna, so
      #   we can share it across applications?
      return FILEPATH_MATCH.match(path)[:bucket_name]
    end

    def self.extract_file_path_from_path(path)
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]
      # Should this be in a central gem, like etna, so
      #   we can share it across applications?
      return FILEPATH_MATCH.match(path)[:file_path]
    end

    def validate_access_to_buckets(user_authorized_bucket_names)
      raise Etna::Forbidden, "Cannot access the source bucket #{source_bucket_name}" unless user_authorized_bucket_names.include? source_bucket_name
      raise Etna::Forbidden, "Cannot access the destination bucket #{dest_bucket_name}" unless user_authorized_bucket_names.include? dest_bucket_name
    end

    def source_bucket_name
      Metis::Revision.extract_bucket_name_from_path(@source)
    end

    def dest_bucket_name
      Metis::Revision.extract_bucket_name_from_path(@dest)
    end

    def source_file_path
      Metis::Revision.extract_file_path_from_path(@source)
    end

    def dest_file_path
      Metis::Revision.extract_file_path_from_path(@dest)
    end

    def source_file
      Metis::File.from_path(
        get_bucket(@source),
        source_file_path)
    end

    private

    def valid_file_path?(path)
      FILEPATH_MATCH.match(path)
    end

    def get_bucket(path)
      Metis::Bucket.where(
        project_name: FILEPATH_MATCH.match(path)[:project_name],
        name: FILEPATH_MATCH.match(path)[:bucket_name]
      ).first
    end

    def self.path_from_parts(project_name, bucket_name, file_path)
      "metis://#{project_name}/#{bucket_name}/#{file_path}"
    end
  end
end
