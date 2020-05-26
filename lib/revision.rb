require_relative 'path'

class Metis
  class Revision
    attr_reader :source, :dest
    def initialize(params)
      @source = Metis::Path.new(params[:source])
      @dest = Metis::Path.new(params[:dest])
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

    def valid? (validation_type=nil, user_authorized_bucket_names=nil)
      return false unless validation_type
      case validation_type
      when 'source_path'
        return @source.valid?
      when 'source_bucket_access'
        return false unless user_authorized_bucket_names
        return user_authorized_bucket_names.include? @source.bucket_name
      end
      return false
    end
  end
end
