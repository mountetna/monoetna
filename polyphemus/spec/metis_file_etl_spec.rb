describe Polyphemus::MetisFileEtl do
  class TestEtl < Polyphemus::MetisFileEtl
    def initialize(**args)
      super(project_bucket_pairs: [['ipi', 'data']], limit: 2, **args)
    end

    def prepare_find_request(cursor, find_request)
      super

      find_request.add_param(Etna::Clients::Metis::FindParam.new(
          type: 'file',
          attribute: 'name',
          predicate: 'glob',
          value: '*',
      ))
    end

    attr_reader :process_calls
    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  # To re-record this test, delete the metis_file_etl2.e2e.yml cassette and TOKEN=xxxxx to your environment before running.
  it 'should process metis files, and support reset' do
    VCR.use_cassette('metis_file_etl.e2e') do
      metis_client = Etna::Clients::Metis.new(host: 'https://metis.development.local', token: ENV['TOKEN'] || 'test-token', persistent: false)
      etl = TestEtl.new(metis_client: metis_client)
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
    all_files = []
    etl.process_calls.each do |call|
      files = call.last
      all_files.push(*files)
    end

    # Do not duplicate any yields
    expect(all_files).to eq(all_files.uniq { |f| f.file_path })

    # Ensure all data is eventually yielded
    expect(all_files.map(&:file_path)).to eq(etl.metis_client.find(Etna::Clients::Metis::FindRequest.new(
        project_name: 'ipi',
        bucket_name: 'data',
        params: [
            Etna::Clients::Metis::FindParam.new(
                type: 'file',
                attribute: 'name',
                predicate: 'glob',
                value: '*',
            )
        ]
    )).files.all.map(&:file_path))
  end
end