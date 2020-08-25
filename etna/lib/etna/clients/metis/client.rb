require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Metis
      def initialize(host:, token:)
        raise 'Metis client configuration is missing host.' unless host
        raise 'Metis client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(host, token)
      end

      def list_all_folders(list_all_folders_request)
        FoldersResponse.new(
          @etna_client.folder_list_all_folders(list_all_folders_request.to_h))
      end

      def list_folder(list_folder_request)
        FoldersAndFilesResponse.new(
          @etna_client.folder_list(list_folder_request.to_h))
      end

      def rename_folder(rename_folder_request)
        FoldersResponse.new(
          @etna_client.folder_rename(rename_folder_request.to_h))
      end

      def create_folder(create_folder_request)
        FoldersResponse.new(
          @etna_client.folder_create(create_folder_request.to_h))
      end

      # Convenience methods ... may belong in a different class?
      def fetch_folders(project_name, bucket_name)
        list_all_folders(
          Etna::Clients::Metis::ListFoldersRequest.new(
            project_name: project_name, bucket_name: bucket_name)).folders
      end

      def parent_folder_path(folder_path)
        folder_path.split('/')[0..-2].join('/')
      end

      def create_parent_folder(project_name, bucket_name, folder_path)
        parent_path = parent_folder_path(folder_path)

        create_folder(
          Etna::Clients::Metis::CreateFolderRequest.new(
            project_name: project_name,
            bucket_name: bucket_name,
            folder_path: parent_path
        ))
      end

      def parent_exists?(project_name, bucket_name, folder_path)
        # NOTE: this doesn't test if the folder_path itself exists
        #   This can be confusing for root folders, because
        #       they have no parents, so you don't need
        #       to create anything.
        parent_path = parent_folder_path(folder_path)
        return true if parent_path.empty? # root folder

        # returns 422 if the folder_path does not exist
        begin
            list_folder(
              Etna::Clients::Metis::ListFolderRequest.new(
                project_name: project_name,
                bucket_name: bucket_name,
                folder_path: parent_path
            ))
        rescue Etna::Error => e
            return false if e.status == 422
            raise
        end
        return true
      end

      def move_folder(project_name, source_bucket_name, source_folder_path, dest_bucket_name, dest_folder_path)
        create_parent_folder(project_name, dest_bucket_name, dest_folder_path) if !parent_exists?(project_name, dest_bucket_name, dest_folder_path)

        rename_folder(
          Etna::Clients::Metis::RenameFolderRequest.new(
            bucket_name: source_bucket_name,
            project_name: project_name,
            folder_path: source_folder_path,
            new_bucket_name: dest_bucket_name,
            new_folder_path: dest_folder_path
        ))
      end
    end
  end
end
