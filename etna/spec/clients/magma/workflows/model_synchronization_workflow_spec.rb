describe Etna::Clients::Magma::ModelSynchronizationWorkflow do
  let(:workflow) do
    Etna::Clients::Magma::ModelSynchronizationWorkflow.from_api_source(source_project: 'mvir1', source_client: source_client)
  end
end