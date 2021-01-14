require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Metis
      class SyncMetisDataWorkflow < Struct.new(:metis_client, :filesystem, :project_name, :bucket_name,
          :logger, :skip_tmpdir, keyword_init: true)
        def copy_directory(src, dest, root = dest, tmpdir = nil)
          own_tmpdir = tmpdir.nil? && !skip_tmpdir
          if own_tmpdir
            tmpdir = filesystem.tmpdir
          end

          begin
            response = metis_client.list_folder(ListFolderRequest.new(project_name: project_name, bucket_name: bucket_name, folder_path: src))

            response.files.all.each do |file|
              logger&.info("Copying file #{file.file_path} (#{Etna::Formatting.as_size(file.size)})")
              copy_file(bin_root_dir: root, tmpdir: tmpdir, dest: ::File.join(dest, file.file_name), url: file.download_url)
            end

            response.folders.all.each do |folder|
              copy_directory(::File.join(src, folder.folder_name), ::File.join(dest, folder.folder_name), root, tmpdir)
            end
          ensure
            filesystem.rm_rf(tmpdir) if own_tmpdir
          end
        end

        def bin_file_name(etag:)
          "bin/#{etag}"
        end

        def with_maybe_intermediate_tmp_dest(bin_file_name:, tmpdir:, dest_file_name:, &block)
          filesystem.mkdir_p(::File.dirname(dest_file_name))
          if tmpdir.nil?
            yield dest_file_name
          else
            tmp_file = ::File.join(tmpdir, ::File.basename(bin_file_name))
            yield tmp_file
            filesystem.mv(tmp_file, dest_file_name)
          end
        end

        def copy_file(bin_root_dir:, tmpdir:, dest:, url:, stub: false)
          metadata = metis_client.file_metadata(url)
          etag = metadata[:etag]
          size = metadata[:size]

          dest_bin_file = ::File.join(bin_root_dir, bin_file_name(etag: etag))
          # Already materialized, continue
          if filesystem.exist?(dest_bin_file)
            return
          end

          with_maybe_intermediate_tmp_dest(bin_file_name: dest_bin_file, tmpdir: tmpdir, dest_file_name: dest) do |tmp_file|
            upload_timings = []
            upload_amount = 0
            last_rate = 0.00001

            filesystem.with_writeable(tmp_file, "w", size_hint: size) do |io|
              if stub
                io.write("(stub) #{size} bytes")
              else
                metis_client.download_file(url) do |chunk|
                  io.write(chunk)

                  upload_timings << [chunk.length, Time.now.to_f]
                  upload_amount += chunk.length

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
                    logger&.info("Uploading #{Etna::Formatting.as_size(rate)} per second")

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

          filesystem.mkdir_p(::File.dirname(dest_bin_file))
          filesystem.with_writeable(dest_bin_file, 'w', size_hint: 0) do |io|
            # empty file marking that this etag has been moved, to save a future write.
          end
        end
      end
    end
  end
end
