require 'ostruct'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Metis
      class SyncMetisDataWorkflow < Struct.new(:metis_client, :filesystem, :project_name, :bucket_name, :logger, keyword_init: true)
        DOWNLOAD_REGEX = /^https:\/\/[^\/]*\/(?<project_name>.*)\/download\/(?<bucket_name>.*)\/(?<file_path>[^\?]*).*$/

        def copy_directory(src, dest, root = dest)
          response = metis_client.list_folder(ListFolderRequest.new(project_name: project_name, bucket_name: bucket_name, folder_path: src))

          response.files.all.each do |file|
            copy_file(dest: ::File.join(dest, file.file_name), url: file.download_url)
          end

          response.folders.all.each do |folder|
            copy_directory(::File.join(src, folder.folder_name), ::File.join(dest, folder.folder_name), root)
          end
        end

        def copy_file(dest:, url:, stub: false)
          # This does not work due to the magma bucket's restrictions, but if it did work, it'd be super sweet.
          # url_match = DOWNLOAD_REGEX.match(url)
          #
          # if filesystem.instance_of?(Etna::Filesystem::Metis) && !url_match.nil?
          #   bucket_name = url_match[:bucket_name]
          #   project_name = url_match[:project_name]
          #   file_path = url_match[:file_path]
          #
          #   metis_client.copy_files(
          #     Etna::Clients::Metis::CopyFilesRequest.new(
          #       project_name: project_name,
          #       revisions: [
          #         Etna::Clients::Metis::CopyRevision.new(
          #           source: "metis://#{project_name}/#{bucket_name}/#{file_path}",
          #           dest: "metis://#{filesystem.project_name}/#{filesystem.bucket_name}#{dest}",
          #         )
          #       ]
          #     )
          #   )
          # end

          metadata = metis_client.file_metadata(url)
          size = metadata[:size]
          etag = metadata[:etag]

          logger&.info("Copying file #{url} (#{Etna::Formatting.as_size(size)})")

          filesystem.with_writeable(dest, "w", size_hint: size) do |io|
            if stub
              io.write("(stub) #{size} bytes")
            else
              metis_client.download_file(url) do |chunk|
                io.write(chunk)
              end
            end
          end
        end
      end
    end
  end
end
