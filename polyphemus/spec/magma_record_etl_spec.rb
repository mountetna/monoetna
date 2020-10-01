describe Polyphemus::MagmaRecordEtl do
  class TestMagmaEtl < Polyphemus::MagmaRecordEtl
    def initialize(**args)
      super(project_model_pairs: [['mvir1', 'timepoint']], limit: 2, **args)
    end

    def prepare_retrieve_request(cursor, retrieve_request)
      super
      # Only fetch the identifiers of the records being scanned.
      retrieve_request.attribute_names = [ 'name', 'identifier', 'updated_at' ]
    end

    attr_reader :process_calls
    def process(cursor, records)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, records.dup]
    end
  end

  # To re-record this test, delete the magma_record_etl.e2e.yml cassette and TOKEN=xxxxx to your environment before running.
  it 'should process magma records, and support reset' do
    VCR.use_cassette('magma_record_etl.e2e') do
      magma_client = Etna::Clients::Magma.new(host: 'https://magma.development.local', token: ENV['TOKEN'] || 'test-token', persistent: false)
      etl = TestMagmaEtl.new(magma_client: magma_client)
      etl.execute('run')
      expect(etl.process_calls.length).to_not eq(0)
      expect(etl.process_calls.length).to_not eq(1)
      check_process_calls(etl)

      etl.process_calls.clear
      etl.execute('run')
      expect(etl.process_calls.length).to eq(0)

      etl.execute('reset')

      etl.process_calls.clear
      etl.execute('run')
      expect(etl.process_calls.length).to_not eq(0)
      check_process_calls(etl)
    end
  end

  def check_process_calls(etl)
    all_processed_records = []
    etl.process_calls.each do |call|
      records = call.last
      all_processed_records.push(*records)
    end

    # Do not duplicate any yields
    expect(all_processed_records).to eq(all_processed_records.uniq)

    # Ensure all data is eventually yielded
    remote_records_ids = etl.magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
        project_name: 'mvir1',
        model_name: 'timepoint',
        record_names: 'all',
    )).models.model('timepoint').documents.document_keys.sort

    expect(all_processed_records.map { |r| r.keys.first }.sort).to eq(remote_records_ids)
  end
end