# Wrap the query_builder and figure out from params
#    what the injected @base_query should be.

require_relative './query_builder'

class Metis
  class Query
    def initialize(project_name:, bucket:, params:)
      @params = params

      @files_base_query = Metis::File.where(
        project_name: @params[:project_name],
        bucket: bucket
      )
      @folders_base_query = Metis::Folder.where(
        project_name: @params[:project_name],
        bucket: bucket
      )

      files_params = filter_params('file')
      folders_params = filter_params('folder')

      @files_query = Metis::QueryBuilder.new(
        @files_base_query,
        files_params
      )

      @folders_query = Metis::QueryBuilder.new(
        @folders_base_query,
        folders_params
      )
    end

    def execute
      {
        files: @files_query.build.all,
        folders: @folders_query.build.all
      }
    end

    private

    def filter_params(type)
      @params.filter { |param|
        !param.key?(:type) || param[:type] == type
      }
    end
  end
end