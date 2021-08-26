require_relative "../../../ipi/ipi_helper"
require_relative "../../../helpers"

class Polyphemus::IpiRnaSeqProcessorBase
  include WithEtnaClients
  include WithLogger

  def initialize
    @helper = IpiHelper.new
  end

  def process(cursor, files)
    raise "Subclasses should implement this method"
  end

  def download_files(files, &block)
    files.each do |attribute_file|
      Tempfile.create do |tmp|
        metis_client.download_file(attribute_file) do |chunk|
          tmp << chunk
        end

        tmp.rewind

        yield [attribute_file, tmp]
      end
    end
  end

  def update_for_cursor(cursor, &block)
    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: cursor[:project_name],
    )

    yield update_request

    magma_client.update_json(update_request) if update_request.revisions.keys.length > 0
  end
end
