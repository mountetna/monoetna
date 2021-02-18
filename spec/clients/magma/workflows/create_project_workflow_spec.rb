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
        tok = 'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6InNhdXJhYmguYXN0aGFuYUB1Y3NmLmVkdSIsImZpcnN0IjoiU2F1cmFiaCIsImxhc3QiOiJBc3RoYW5hIiwicGVybSI6IkE6YWRtaW5pc3RyYXRpb247YTpkc2NvbGFiLGR1bmxhcCxpcGkiLCJleHAiOjE2MTM3MDQ0MDZ9.aI3Tj79lQXFDghrfj7lcCr9FUPKC0Em80KUJsHZ__92HxcqVHJbVG1M1uOSDAYNvz8gn21g29atPDC1kRY2WX39tTlVvmrTG82TCxk_kEcsYPYvcvCrWcz8SHu8tHFwJdM1A19YOuAsed6Do5RaupUNikqwDYld1ZjxtyyEhFNVvClk5EtLKHPbfS4XLieJNx0ab6Eob-XTk_ADRSuczeQb4klBVDZ6OA7ynzeJLd72DHF4RArGcrVWM0bfVOMEbbDjhGEzK324N29CFsoCBOylUkF4SiP92_Eb5h1Kgq2pbS3vtyNWGXy-p6gPNV2rFoK8mWhO3sSsnljqO0MWMYg'

        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: tok,
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
