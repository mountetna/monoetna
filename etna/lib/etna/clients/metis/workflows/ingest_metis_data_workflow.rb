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
        def copy_files(srcs)
          srcs.each do |src|
            next unless ingest_filesystem.exist?(src)

            logger&.info("Copying file #{src} (#{Etna::Formatting.as_size(ingest_filesystem.stat(src).size)})")

            # For ingestion triage, just copy over the exact path + filename.
            copy_file(dest: src, src: src)
          end
        end

        def copy_file(dest:, src:, stub: false)
          ingest_filesystem.with_readable(src, "r") do |file|
            metis_filesystem.do_streaming_upload(file, dest, file.size)
          end
        end
      end
    end
  end
end
