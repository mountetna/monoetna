require 'webmock/rspec'
require 'json'
require_relative '../../../../lib/etna/clients/janus'

describe Etna::Clients::Magma::CreateProjectWorkflow do
  describe "e2e" do
    it 'can create a new project in janus and magma' do
      configure_etna_yml

      allow(STDIN).to receive(:gets).and_return('n')

      Timecop.freeze('2000-01-01') do
        VCR.use_cassette('create_project_workflow.e2e') do
          # Change this when re-running a cassette to ensure a new project is being created in your environment
          test_project = "test_create_project_aai"

          # This is an expired development token and is safe to make public, does not leak anything about production or staging values
          # and cannot be used in a sensitive way.
          # This test depends on breaking down a tok to get user info, so it's important we actually set a token
          tok = 'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsIm5hbWUiOiJEZXZlbG9wZXIgTGFzdE5hbWUiLCJwZXJtIjoiQTp0ZXN0LXByb2plY3Q7YTphZG1pbmlzdHJhdGlvbiIsImV4cCI6ODY0MDE2MDgxMzY2MzV9.cL1PAZHf1IfegxkdXDUKSzhdgO5c_wM_veZsp7pyg28uxyXGs9CKWRNyh2zAIrdHF40Hb6rLDj0mr1itaUtILcqTVEss2SB5xhHly3gUdKAkjVACB6FNltUrBhobBlhKxbL_-PS0SI5eRrIHQdQ1iO5wMAcjHGqSh3zYC7SsNHiW9WjwvRlJ2SJSOzq_kJB5G_wHI1ZqsSHuMarc22ztQ_lU9CrqrX4csGkNm2bouYJginXggCVgI34Cwj6nmBWHte0u8YOBUxk2OO1Nln4YFaMFAlPvnD9KzU2CQAO93IDt5bAXWURWeYNYA0A6kcx-GLmvx1RbcnSZR5z_PoCmiQ'

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
end
