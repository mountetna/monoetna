describe Etna::Clients::Magma::AddProjectModelsWorkflow do
  describe "e2e" do
    it 'can sync a full project from point A to point B' do
      configure_etna_yml

      VCR.use_cassette('add_project_models_workflow-full-project.e2e') do
        # Change this name when re-recording the cassette file to ensure a new project is synced
        test_project = "test_add_project_models_workflow_full_aad"

        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: ENV['TOKEN'] || 'token', persistent: false, ignore_ssl: true,
        )

        magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: test_project,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))

        # Ensure this is a brand new project we are testing with.
        expect(magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: test_project)).models.model_keys).to eql(['project'])

        workflow = Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)

        io = StringIO.new
        workflow.write_models_templats_csv(io, 'mvir1')

        # io.rewind
        # puts io.read

        io.rewind
        models = workflow.prepare_models_from_csv(io) do |err|
          raise err
        end

        workflow.synchronize_to_server(models, test_project) do |update|
          pp update
        end
      end
    end
  end
end
