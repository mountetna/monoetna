require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Metis
      class MetisDownloadWorkflow < Struct.new(:metis_client, :project_name, :bucket_name, :max_attempts, keyword_init: true)

        def initialize(args)
          super({max_attempts: 3}.update(args))
        end

        # TODO: Might be possible to use range headers to select and resume downloads on failure in the future.
        def do_download(dest_file_or_io, metis_file, &block)
          size = metis_file.size
          completed = 0.0
          start = Time.now

          unless dest_file_or_io.is_a?(IO) || dest_file_or_io.is_a?(StringIO)
            ::File.open(dest_file_or_io, 'w') do |io|
              return do_download(dest_file_or_io, metis_file, &block)
            end
          end

          metis_client.download_file(metis_file) do |chunk|
            dest_file_or_io.write chunk
            completed += chunk.size

            block.call([
                :progress,
                size == 0 ? 1 : completed / size,
                (completed / (Time.now - start)).round(2),
            ]) unless block.nil?
          end
        end
      end
    end
  end
end
