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
            logger&.info("Copying file #{file.file_path} (#{Etna::Formatting.as_size(file.size)})")
            copy_file(dest: ::File.join(dest, file.file_name), url: file.download_url)
          end

          response.folders.all.each do |folder|
            copy_directory(::File.join(src, folder.folder_name), ::File.join(dest, folder.folder_name), root)
          end
        end

        def copy_file(dest:, url:, stub: false)
          url_match = DOWNLOAD_REGEX.match(url)

          if filesystem.instance_of?(Etna::Filesystem::Metis) && !url_match.nil?
            bucket_name = url_match[:bucket_name]
            project_name = url_match[:project_name]
            file_path = url_match[:file_path]

            metis_client.copy_files(
              Etna::Clients::Metis::CopyFilesRequest.new(
                project_name: project_name,
                revisions: [
                  Etna::Clients::Metis::CopyRevision.new(
                    source: "metis://#{project_name}/#{bucket_name}/#{file_path}",
                    dest: "metis://#{filesystem.project_name}/#{filesystem.bucket_name}/#{dest}",
                  )
                ]
              )
            )

            return
          end

          metadata = metis_client.file_metadata(url)
          size = metadata[:size]

          tmp_file = dest
          upload_timings = []
          upload_amount = 0
          last_rate = 0.00001
          remaining = size

          filesystem.with_writeable(tmp_file, "w", size_hint: size) do |io|
            if stub
              io.write("(stub) #{size} bytes")
            else
              metis_client.download_file(url) do |chunk|
                io.write(chunk)

                upload_timings << [chunk.length, Time.now.to_f]
                upload_amount += chunk.length
                remaining -= chunk.length

                if upload_timings.length > 150
                  s, _ = upload_timings.shift
                  upload_amount -= s
                end

                _, start_time = upload_timings.first
                _, end_time = upload_timings.last

                if start_time == end_time
                  next
                end

                rate = upload_amount / (end_time - start_time)

                if rate / last_rate > 1.3 || rate / last_rate < 0.7
                  logger&.info("Uploading #{Etna::Formatting.as_size(rate)} per second, #{Etna::Formatting.as_size(remaining)} remaining")

                  if rate == 0
                    last_rate = 0.0001
                  else
                    last_rate = rate
                  end
                end
              end

            end
          end
        end
      end
    end
  end
end
