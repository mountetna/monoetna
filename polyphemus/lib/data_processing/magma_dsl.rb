module MagmaDsl
  def load_magma_record(table, identifier)
    workflow = Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: @magma_client, project_name: @project_name)
    workflow.lookup_record(table, identifier)
  end
end