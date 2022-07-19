describe Polyphemus::MagmaRecordScriptEtl::Scripts do
  describe Polyphemus::MagmaRecordScriptEtl::Scripts::IpiPatientFlowjo do
    xit 'runs without exception' do
      VCR.use_cassette('ipi+patient+flowjo.e2e') do
        magma_client = Etna::Clients::Magma.new(host: 'https://magma.ucsf.edu', token: ENV['TOKEN'] || TEST_TOKEN)
        metis_client = Etna::Clients::Metis.new(host: 'https://metis.ucsf.edu', token: ENV['TOKEN'] || TEST_TOKEN)
        crud = Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: 'ipi')

        etl = described_class.new
        allow(etl).to receive(:magma_client).and_return(magma_client)
        allow(etl).to receive(:metis_client).and_return(metis_client)

        etl.process(nil, [crud.lookup_record('patient', 'IPIGSTR001')])
      end
    end
  end
end
