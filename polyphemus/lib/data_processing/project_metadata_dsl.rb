module ProjectMetadataDsl
  def models
    @models ||= begin
      magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
    end
  end

  def record_identifier(record = @record, model_name = self.model_name)
    record[models.model(model_name).template.identifier]
  end
end