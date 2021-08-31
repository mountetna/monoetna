# Returns a set of Metis files, if they are in a watch folder
require_relative "./metis_file_etl"

class Polyphemus
  class MetisFileInWatchFolderEtl < MetisFileEtl
    def initialize(project_bucket_pairs:, metis_client: nil, limit: 20, file_name_params: {}, watch_type:)
      super(
        project_bucket_pairs: project_bucket_pairs,
        metis_client: metis_client,
        limit: limit,
        file_name_params: file_name_params,
      )
      @watch_type = watch_type
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
    end

    def watch_folder_ids(cursor)
      @watch_folder_ids ||= Polyphemus::WatchFolder.where(
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        watch_type: @watch_type,
      ).all.map do |folder|
        folder.folder_id
      end
    end

    # def execute_request(find_request, i)
    #   # Only pass along files that are in a watch folder
    #   # Because this filters out files, we also need to iterate here
    #   #   and increase the request limit until we get actual
    #   #   results, otherwise the ETL could get into a state
    #   #   where it never does work, because the first @limit
    #   #   files from Metis aren't in watch folders.
    #   files_in_watch_folders = []

    #   loop do
    #     files = super(find_request, i)

    #     break if files.empty?

    #     files_in_watch_folders = files.select do |file|
    #       watch_folder_paths(
    #         find_request.project_name, find_request.bucket_name
    #       ).include?(::File.dirname(file.file_path))
    #     end

    #     break unless files_in_watch_folders.empty?
    #     i += 1
    #   end

    #   files_in_watch_folders
    # end

    # def watch_folder_paths(project_name, bucket_name)
    #   @watch_folder_paths ||= Polyphemus::WatchFolder.where(
    #     project_name: project_name,
    #     bucket_name: bucket_name,
    #     watch_type: @watch_type,
    #   ).all.map do |folder|
    #     folder.folder_path
    #   end
    # end
  end
end
