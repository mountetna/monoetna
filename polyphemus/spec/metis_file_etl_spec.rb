describe Polyphemus::MetisFileEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestMetisEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands['test_metis_etl']
  end

  def run_etl_command(*args)
    cmd, args, kwds = etl_command.find_command(*args)
    cmd.execute(*args, **kwds)
  end

  def setup_client(metis_client)
    allow(etl_executor.subcommands['run'].etl).to receive(:metis_client).and_return(metis_client)
    allow(etl_executor.subcommands['reset'].etl).to receive(:metis_client).and_return(metis_client)
  end

  class TestMetisEtl < Polyphemus::MetisFileEtl
    def initialize(**args)
      super(project_bucket_pairs: [['ipi', 'data']], limit: 2, **args)
    end

    def prepare_find_request(cursor, find_request)
      super
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
      metis_client = Etna::Clients::Metis.new(host: 'https://metis.development.local', token: ENV['TOKEN'] || TEST_TOKEN)
      setup_client(metis_client)

      etl = etl_executor.subcommands['run'].etl

      run_etl_command('test_metis_etl', 'run')
      expect(etl.process_calls.length).to_not eq(0)
      expect(etl.process_calls.length).to_not eq(1)
      check_process_calls(etl)

      etl.process_calls.clear
      run_etl_command('test_metis_etl', 'run')
      expect(etl.process_calls.length).to eq(0)

      run_etl_command('test_metis_etl', 'reset')

      etl.process_calls.clear
      run_etl_command('test_metis_etl', 'run')
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
                predicate: '=~',
                value: '%',
            )
        ]
    )).files.all.map(&:file_path))
  end
end