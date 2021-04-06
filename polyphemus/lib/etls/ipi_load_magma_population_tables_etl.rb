require_relative "../magma_record_etl"

class Polyphemus::IpiLoadMagmaPopulationTablesEtl < Polyphemus::MagmaRecordEtl
  def initialize
    super(project_model_pairs: [["ipi", "patient"]], attribute_names: ["ipi_number", "updated_at"])
  end

  def process(cursor, records)
    project_name = "ipi"
    record_names = records.map { |r| r.keys.first }
    logger.info("Processing population tables for patients #{record_names.join(", ")}...")

    ipi_flowjo_script = File.join(File.dirname(__FILE__), "magma", "scripts", "ipi+patient+flowjo.rb")

    runner = Polyphemus::MagmaEtlScriptRunner.new(ipi_flowjo_script)

    crud = Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)

    processed_wsp_updated_recently?(cursor, record_names).each do |record_id|
      record = crud.lookup_record(runner.model_name, record_id)
      if record.nil?
        raise "Could not find #{runner.model_name} by id #{record_id} in project #{project_name}"
      end

      begin
        runner.run(record, magma_client: magma_client, commit: true)
      rescue => e
        logger.error("#{e.message}\n#{e.backtrace}")
      end
    end

    logger.info("Done")
  end

  def processed_wsp_updated_recently?(cursor, record_names)
    magma_client.query(Etna::Clients::Magma::QueryRequest.new(
      project_name: "ipi",
      query: ["patient", ["ipi_number", "::in", record_names],
              "::all", "flojo_file_processed", "::updated_at"],
    )).answer.map do |result|
      record_name = result[0]
      record_file_updated_at = result[1]
      cursor_record = cursor[:seen_ids].find { |s|
        s[0] == record_name
      }

      # Only run the population loader for the patient
      #   if they weren't in the cursor's
      #   last batch of seen_ids, or if the file has been updated since
      #   the cursor's last recorded updated_at for the given patient.
      if !record_file_updated_at
        nil
      elsif !cursor_record
        record_name
      else
        # Note that in Magma, the patient record is updated
        #   slightly after the file is updated on Metis,
        #   so we're conservative and compare against
        #   the min of the cursor's updated_at or the seen_id record.
        # This lets updated files be identified both for the initial
        #   run of this ETL as well as subsequent runs.
        last_updated_at = [Time.at(cursor_record[1]), cursor.updated_at].min
        record_name if Time.parse(record_file_updated_at) >= last_updated_at
      end
    end.compact
  end
end
