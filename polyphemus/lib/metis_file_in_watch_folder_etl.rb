# Returns a set of Metis files, if they are in a watch folder
require_relative "./metis_file_etl"

class Polyphemus
  class MetisFileInWatchFolderEtl < MetisFileEtl
    def initialize(project_bucket_pairs:, metis_client: nil, limit: 20, watch_type:, hide_paths: true)
      super(
        project_bucket_pairs: project_bucket_pairs,
        metis_client: metis_client,
        limit: limit,
        file_name_params: {},
      )
      @watch_type = watch_type
      @hide_paths = hide_paths
    end

    def prepare_find_request(cursor, find_request)
      super

      # Here we update the find request with the folder_ids from
      #   the watch folders
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
        type: "file",
        attribute: "folder_id",
        predicate: "in",
        value: watch_folder_ids(cursor),
      ))
      find_request.hide_paths = @hide_paths
    end

    def watch_folder_ids(cursor)
      @watch_folder_ids ||= Polyphemus::WatchFolder.where(
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        watch_type: @watch_type,
      ).all.map do |folder|
        folder.metis_id
      end
    end
  end
end
