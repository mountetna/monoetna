require_relative "../magma_record_etl"

class Polyphemus::IpiLoadMagmaPopulationTablesEtl < Polyphemus::MagmaRecordEtl
  def initialize
    super(project_model_pairs: [["ipi", "patient"]], attribute_names: ["ipi_number", "updated_at"])
  end

  def process(cursor, records)
    project_name = "ipi"
    record_names = records.map { |r| r.keys.first }
    logger.info("Processing patients #{record_names.join(", ")}...")

    ipi_flowjo_script = File.join(File.dirname(__FILE__), "magma", "scripts", "ipi+patient+flowjo.rb")

    runner = Polyphemus::MagmaEtlScriptRunner.new(ipi_flowjo_script)

    crud = Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)

    record_names.each do |record_id|
      record = crud.lookup_record(runner.model_name, record_id)
      if record.nil?
        raise "Could not find #{runner.model_name} by id #{record_id} in project #{project_name}"
      end

      runner.run(record, magma_client: magma_client, commit: true)
    end

    logger.info("Done")
  end
end
