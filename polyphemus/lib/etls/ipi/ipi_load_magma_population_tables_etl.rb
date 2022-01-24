require_relative "../../magma_record_file_etl"

class Polyphemus::IpiLoadMagmaPopulationTablesEtl < Polyphemus::MagmaRecordFileEtl
  def initialize(cursor_env: {}, scanner: nil)
    super(
      project_model_pairs: [["ipi", "patient"]],
      attribute_names: ["updated_at", "flojo_file_processed"],
      cursor_env: cursor_env,
      scanner: scanner
    )
  end

  def process(cursor, records)
    project_name = "ipi"
    record_names = records.map { |r| r.keys.first }
    logger.info("Processing population tables for patients #{record_names.join(", ")}...")

    ipi_flowjo_script = File.join(File.dirname(__FILE__), "magma", "scripts", "ipi+patient+flowjo.rb")

    runner = Polyphemus::MagmaEtlScriptRunner.new(ipi_flowjo_script)

    crud = Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)

    record_names.each do |record_id|
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
end
