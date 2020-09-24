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
        def do_download(dest_file, metis_file, &block)
          size = metis_file.size
          completed = 0.0
          start = Time.now

          ::File.open(dest_file, "w") do |io|
            metis_client.download_file(metis_file) do |chunk|
              io.write chunk
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
end
