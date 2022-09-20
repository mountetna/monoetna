require_relative '../magma_record_etl'

class Polyphemus::MaterializeGneMagmaRecordsEtl < Polyphemus::MagmaRecordEtl
  def initialize(cursor_env: {}, scanner: build_scanner)
    super(project_model_pairs: [['mvir1', 'patient']], attribute_names: ['name', 'updated_at'], cursor_env: cursor_env, scanner: scanner)
  end

  def process(cursor, records)
    record_names = records.map { |r| r.keys.first }
    logger.info("Processing patients #{record_names.join(', ')}...")
    workflow = Etna::Clients::Magma::MaterializeDataWorkflow.new(
        model_attributes_mask: model_attribute_pairs,
        record_names: record_names,
        model_filters: model_filters,
        metis_client: metis_client, magma_client: magma_client, logger: logger,
        project_name: 'mvir1', model_name: 'patient', filesystem: filesystem)

    workflow.materialize_all("Upload")
    logger.info("Done")
  end

  def filesystem
    @filesystem ||= Etna::Filesystem::Metis.new(metis_client: metis_client,
        project_name: 'mvir1', bucket_name: 'GNE_composite')
  end

  private

  def model_filters
    default_model_filters.update(model_filters_override)
  end

  def default_model_filters
    {
        'patient' => 'comet_plus=true'
    }
  end

  def model_filters_override
    Polyphemus.instance.config(:gne_model_filters) || {}
  end

  def model_attribute_pairs
    default_model_attribute_pairs.update(model_attributes_override)
  end

  def default_model_attribute_pairs
    result = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
        project_name: 'mvir1',
        model_name: 'all',
        attribute_names: 'identifier',
        record_names: 'all'
    )).models.model_keys.map do |model_name|
        [model_name, []]
    end.to_h

    # This ensures the default model_filters can be applied
    result['patient'] << 'comet_plus'

    result
  end

  def model_attributes_override
    Polyphemus.instance.config(:gne_model_attributes) || {}
  end
end