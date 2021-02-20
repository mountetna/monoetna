require_relative '../metis_file_etl'
require 'concurrent'

class Polyphemus::SyncGneMetisFilesEtl < Polyphemus::MetisFileEtl
  BUCKET = 'GNE_redacted_data'

  def initialize
    super(project_bucket_pairs: [['mvir1', BUCKET]])
  end

  def process(cursor, files)
    workflow = Etna::Clients::Metis::SyncMetisDataWorkflow.new(
        metis_client: metis_client, logger: logger,
        project_name: 'mvir1', bucket_name: BUCKET,
        filesystem: filesystem)

    concurrency = 10
    root = "Upload/processed/#{BUCKET}"

    semaphore = Concurrent::Semaphore.new(concurrency)
    errors = Queue.new

    files.each do |file|
      begin
        if (error = errors.pop(true))
          raise error
        end
      rescue ThreadError
      end

      semaphore.acquire
      Thread.new do
        begin
          file_path = file.file_path
          file_path = file_path.gsub(/^data\//, "")
          workflow.copy_file(dest: ::File.join(root, file_path), url: file.download_url)
        rescue => e
          errors << e
        ensure
          semaphore.release
        end
      end
    end

    semaphore.acquire(concurrency)

    begin
      if (error = errors.pop(true))
        raise error
      end
    rescue ThreadError
    end

    logger.info("Done")
  end

  def filesystem
    @filesystem ||= Etna::Filesystem::Metis.new(metis_client: metis_client,
        project_name: 'mvir1', bucket_name: 'GNE_composite')
  end
  # def filesystem
    # aspera_comet = Polyphemus.instance.config(:aspera_comet)
    # @filesystem ||= Etna::Filesystem::GeneAsperaCliFilesystem.new(**aspera_comet)
  # end
end