require 'webmock/rspec'
require 'json'
require_relative '../../../../lib/etna/clients/janus'

describe Etna::Clients::Magma::CreateProjectWorkflow do
  describe "e2e" do
    it 'can create a new project in janus and magma' do
      configure_etna_yml

      allow(STDIN).to receive(:gets).and_return('n')

      VCR.use_cassette('create_project_workflow.e2e') do
        # Change this when re-running a cassette to ensure a new project is being created in your environment
        test_project = "test_create_project_aaf"

        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: TEST_TOKEN,
            ignore_ssl: true,
        )

        janus_client = Etna::Clients::Janus.new(
            host: 'https://janus.development.local',
            token: TEST_TOKEN,
            ignore_ssl: true,
        )

        workflow = Etna::Clients::Magma::CreateProjectWorkflow.new(
            magma_client: magma_client,
            janus_client: janus_client,
            project_name: test_project,
            project_name_full: 'This is a test ' + test_project
        )

        workflow.create!
      end
    end
  end
end
