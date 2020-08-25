module Etna
  module Clients
    class ProjectClientAdapter
      def initialize(project_name, metis_client)
        @project_name = project_name
        @metis_client = metis_client
      end

      def fetch_folders(bucket_name)
        @metis_client.list_all_folders(
          Etna::Clients::Metis::ListFoldersRequest.new(
            project_name: @project_name, bucket_name: bucket_name)).folders
      end

      def parent_folder_path(folder_path)
        folder_path.split('/')[0..-2].join('/')
      end

      def create_parent_folder(bucket_name, folder_path)
        parent_path = parent_folder_path(folder_path)

        @metis_client.create_folder(
          Etna::Clients::Metis::CreateFolderRequest.new(
            project_name: @project_name,
            bucket_name: bucket_name,
            folder_path: parent_path
        ))
      end

      def parent_exists?(bucket_name, folder_path)
        # NOTE: this doesn't test if the folder_path itself exists
        #   This can be confusing for root folders, because
        #       they have no parents, so you don't need
        #       to create anything.
        parent_path = parent_folder_path(folder_path)
        return true if parent_path.empty? # root folder

        # returns 422 if the folder_path does not exist
        begin
          @metis_client.list_folder(
              Etna::Clients::Metis::ListFolderRequest.new(
                project_name: @project_name,
                bucket_name: bucket_name,
                folder_path: parent_path
            ))
        rescue Etna::Error => e
            return false if e.status == 422
            raise
        end
        return true
      end

      def rename_folder(source_bucket_name, source_folder_path, dest_bucket_name, dest_folder_path)
        create_parent_folder(dest_bucket_name, dest_folder_path) if !parent_exists?(dest_bucket_name, dest_folder_path)

        @metis_client.rename_folder(
          Etna::Clients::Metis::RenameFolderRequest.new(
            bucket_name: source_bucket_name,
            project_name: @project_name,
            folder_path: source_folder_path,
            new_bucket_name: dest_bucket_name,
            new_folder_path: dest_folder_path
        ))
      end
    end
  end
end