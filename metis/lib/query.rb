# Wrap the query_builder and figure out from params
#    what the injected @base_query should be.

require_relative './query_builder'

class Metis
  class Query
    def initialize(project_name:, bucket:, params:)
      @project_name = project_name
      @bucket = bucket
      @params = params

      @files_query = get_query('file')
      @folders_query = get_query('folder')
    end

    def execute
      {
        files: @files_query ? @files_query.build.all : [],
        folders: @folders_query ? @folders_query.build.all : []
      }
    end

    private

    def get_query(query_type)
      base_query = query_type == 'file' ?
        Metis::File.where(
          project_name: @project_name,
          bucket: @bucket
        ) :
        Metis::Folder.where(
          project_name: @project_name,
          bucket: @bucket
        )
      query_params = filter_params(query_type)

      query_params.length > 0 ? Metis::QueryBuilder.new(
        base_query,
        query_params
      ) : nil
    end

    def filter_params(type)
      @params.select { |param|
        !param.key?(:type) || param[:type] == type
      }
    end
  end
end