describe Polyphemus::MetisFileEtl do
  class TestEtl < Polyphemus::MetisFileEtl
    def initialize(**args)
      super(project_bucket_pairs: [['ipi', 'data'], ['ipi', 'waiver_data']], **args)
      @cursor_group.cursors.each { |c| c[:limit] = 2} # Force paging on everything
    end

    attr_reader :process_calls
    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  # To re-record this test, add RERECORD=1 and METIS_TOKEN=xxxxx to your environment before running.
  it 'should run end to end' do
    VCR.use_cassette('metis_file_etl.e2e') do
      metis_client = Etna::Clients::Metis.new(host: 'https://metis.development.local', token: ENV['METIS_TOKEN'])
      etl = TestEtl.new(metis_client: metis_client)
      etl.execute('run')
    end
  end
end