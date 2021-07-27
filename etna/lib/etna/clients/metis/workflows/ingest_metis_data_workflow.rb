require "ostruct"
require "fileutils"
require "tempfile"

module Etna
  module Clients
    class Metis
      class IngestMetisDataWorkflow < Struct.new(:metis_filesystem, :ingest_filesystem, :logger, keyword_init: true)
        # Since we are doing manual triage of files,
        #   do not automatically copy directory trees.
        # srcs must be a list of full paths to files.
        def copy_files(srcs, &block)
          srcs.each do |src|
            if !ingest_filesystem.exist?(src)
              logger&.warn("#{src} does not exist on source filesystem. Skipping.")
              next
            end

            logger&.info("Copying file #{src} (#{Etna::Formatting.as_size(ingest_filesystem.stat(src).size)})")

            # For ingestion triage, just copy over the exact path + filename.
            if block_given?
              yield copy_file(src, src, &block)
            else
              copy_file(src, src)
            end
          end
        end

        def copy_file(dest, src, &block)
          ingest_filesystem.with_readable(src, "r") do |file|
            metis_filesystem.do_streaming_upload(file, dest, file.size)
            yield src if block_given?
          end
        end
      end
    end
  end
end
