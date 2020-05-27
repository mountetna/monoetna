require_relative 'path'

class Metis
  class Revision
    attr_reader :source, :dest, :errors
    def initialize(params)
      @source = Metis::Path.new(params[:source])
      @dest = Metis::Path.new(params[:dest])

      @errors = []
      # TODO: Will have to handle DELETE differently,
      #       when params[:dest] will be nil
    end

    def self.create_from_parts(params)
      Metis::Revision.new({
        source: Metis::Path.path_from_parts(
          params[:source][:project_name],
          params[:source][:bucket_name],
          params[:source][:file_path]),
        dest: Metis::Path.path_from_parts(
          params[:dest][:project_name],
          params[:dest][:bucket_name],
          params[:dest][:file_path])
      })
    end

    def validate (user_authorized_bucket_names)
      @errors.push("Invalid source path: #{@source.path}") unless @source.valid?

      if @source.valid?
        @errors.push("Invalid source bucket: #{@source.bucket_name}") unless @source.bucket

        if @source.bucket
          @errors.push("File #{@source.path} not found") unless @source.file&.has_data?
          @errors.push(
            "Forbidden: no access to source bucket #{@source.bucket.name}"
          ) unless user_authorized_bucket_names.include? @source.bucket.name
        end
      end
    end

    def bucket_names
      # This is kind of weird, but we need the ability to grab
      #   all relevant bucket names, even before validation of
      #   the CopyRevision (in file_controller).
      # So, if @dest doesn't exist, we only return the source
      #   bucket name.
      if @source.valid?
        return [@source.bucket_name]
      end

      return []
    end
  end
end
