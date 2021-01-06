describe Etna::Clients::Magma::ModelSynchronizationWorkflow do
  # Swap this out if you like when you need to re-generate the recording.
  let(:target_project) { "testproject1234111" }
  let(:workflow) do
    @source_client = Etna::Clients::Magma.new(host: "https://magma.ucsf.edu", token: ENV['PROD_TOKEN'] || 'token')
    @target_client = Etna::Clients::Magma.new(host: "https://magma.development.local", token: ENV['LOCAL_TOKEN'] || 'token')

    @target_client.update_model(
        Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: target_project,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))

    Etna::Clients::Magma::ModelSynchronizationWorkflow.from_api_source(
        source_project: 'mvir',
        source_client: @source_client,
        target_client: @target_client,
        target_project: target_project,
    )
  end

  xit 'can synchronize a production model locally' do
    # VCR.use_cassette('model_synchronization_workflow.deep') do
      workflow.ensure_model_tree('project')
      models = @target_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: target_project)).models
    # end
  end
end