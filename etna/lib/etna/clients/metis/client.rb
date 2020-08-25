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
    end
  end
end
