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

        # This is an expired development token and is safe to make public, does not leak anything about production or staging values
        # and cannot be used in a sensitive way.
        # This test depends on breaking down a tok to get user info, so it's important we actually set a token
        tok = 'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsImZpcnN0IjoiRGV2ZWxvcGVyIiwibGFzdCI6IlBlcnNvbiIsInBlcm0iOiJhOmFkbWluaXN0cmF0aW9uLGNvcHJvamVjdF90ZW1wbGF0ZSxpcGksbXZpcjEsdGVzdC1wcm9qZWN0LHRlc3RfY3JlYXRlX3Byb2plY3RfYWFhLHRlc3RfY3JlYXRlX3Byb2plY3RfYWFiLHRlc3RfY3JlYXRlX3Byb2plY3RfYWFjLHRlc3RfbXZpcjEsdGVzdF9tdmlyMix0ZXN0X212aXIzLHRlc3RfbXZpcl9hYmMiLCJleHAiOjg2NDAxNjA4MTM2NjM1fQ.I7eIxpAC7w_KjceR9_QEcWHqQOlwxBAj29s0liV7URNaKSLtTeDLS8eSFa2QCFAqqRzbU3HNDOMopNbKCKaJKFeJC4tWOHaDGYhzN3oQHXmD6vpCB6Ty3fRtz8D3y0sMHNM5ywNopJOPeGaSKUMoLMxvINBXfzFh48IhSs9WEt1myKGwNyvPMZntqZWMTBfXU22FQuVPApUEmhDTYgS5kiPZ2RR0l6yNPmTSwUteStWc4wzJFS73yXI-Lhpbbvza8ih4VU4cfxUso6-V1k3FGFgXMX6O3qEwqTuwLosM8sunn6lyJck_U-DXmpVywU0LcMsi2hq8dLUVTjXEY1qIcw'

        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: tok,
            persistent: false,
            ignore_ssl: true,
        )

        janus_client = Etna::Clients::Janus.new(
            host: 'https://janus.development.local',
            token: tok,
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
