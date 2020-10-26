require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require 'cgi'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Metis
      attr_reader :token
      def initialize(host:, token:, persistent: true, ignore_ssl: false)
        raise 'Metis client configuration is missing host.' unless host
        raise 'Metis client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(
          host,
          token,
          persistent: persistent,
          ignore_ssl: ignore_ssl)

        @token = token
      end

      def list_all_folders(list_all_folders_request = ListFoldersRequest.new)
        FoldersResponse.new(
          @etna_client.folder_list_all_folders(list_all_folders_request.to_h))
      end

      def list_folder(list_folder_request = ListFolderRequest.new)
        FoldersAndFilesResponse.new(
          @etna_client.folder_list(list_folder_request.to_h))
      end

      def ensure_parent_folder_exists(project_name:, bucket_name:, path:)
        create_folder_request = CreateFolderRequest.new(
          project_name: project_name,
          bucket_name: bucket_name,
          folder_path: parent_folder_path(path)
        )
        create_folder(create_folder_request) if !folder_exists?(create_folder_request)
      end

      def rename_folder(rename_folder_request)
        ensure_parent_folder_exists(
          project_name: rename_folder_request.project_name,
          bucket_name: rename_folder_request.new_bucket_name,
          path: rename_folder_request.new_folder_path
        ) if rename_folder_request.create_parent

        FoldersResponse.new(
          @etna_client.folder_rename(rename_folder_request.to_h))
      end

      def create_folder(create_folder_request)
        FoldersResponse.new(
          @etna_client.folder_create(create_folder_request.to_h))
      end

      def find(find_request)
        FoldersAndFilesResponse.new(
          @etna_client.bucket_find(find_request.to_h))
      end

      def download_file(file_or_url = File.new, &block)
        if file_or_url.instance_of?(File)
          download_path =  file_or_url.download_path
        else
          download_path = file_or_url.sub(%r!^https://[^/]*?/!, '/')
        end

        @etna_client.get(download_path) do |response|
          response.read_body(&block)
        end
      end

      def upload_start(upload_start_request = UploadStartRequest.new)
        json = nil
        @etna_client.post(upload_start_request.upload_path, upload_start_request) do |res|
          json = JSON.parse(res.body)
        end

        UploadResponse.new(json)
      end

      def authorize_upload(authorize_upload_request = AuthorizeUploadRequest.new)
        json = nil
        @etna_client.post("/authorize/upload", authorize_upload_request) do |res|
          json = JSON.parse(res.body)
        end

        UploadResponse.new(json)
      end

      def upload_blob(upload_blob_request = UploadBlobRequest.new)
        json = nil
        @etna_client.multipart_post(upload_blob_request.upload_path, upload_blob_request.encode_multipart_content) do |res|
          json = JSON.parse(res.body)
        end

        UploadResponse.new(json)
      end

      def copy_files(copy_files_request)
        FilesResponse.new(
          @etna_client.file_bulk_copy(copy_files_request.to_h))
      end

      def folder_exists?(create_folder_request)
        # NOTE: this doesn't test if the folder_path itself exists
        #   This can be confusing for root folders, because
        #       they have no parents, so you don't need
        #       to create anything.
        return true if create_folder_request.folder_path.empty? # root folder

        # returns 422 if the folder_path does not exist
        begin
          list_folder(
              Etna::Clients::Metis::ListFolderRequest.new(
                project_name: create_folder_request.project_name,
                bucket_name: create_folder_request.bucket_name,
                folder_path: create_folder_request.folder_path
            ))
        rescue Etna::Error => e
            return false if e.status == 422
            raise
        end
        return true
      end

      def folders(project_name:, bucket_name:)
        @folders ||= Hash.new { |h, key|
          h[key] = list_all_folders(
            Etna::Clients::Metis::ListFoldersRequest.new(
              project_name: project_name,
              bucket_name: key
            )).folders.all
        }

        @folders[bucket_name]
      end

      def rename_folders_by_regex(project_name:, source_bucket:, source_folders:, dest_bucket:, regex:)
        found_folders = source_folders.select { |folder|
            folder.folder_path =~ regex
        }

        return if found_folders.length == 0

        found_folders.each { |folder|
            rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
              bucket_name: source_bucket,
              project_name: project_name,
              folder_path: folder.folder_path,
              new_bucket_name: dest_bucket,
              new_folder_path: folder.folder_path,
              create_parent: true)
            )
        }
      end

      private

      def parent_folder_path(folder_path)
        folder_path.split('/')[0..-2].join('/')
      end
    end
  end
end
