require 'etna/clients'
require 'etna/environment_scoped'

module WithEtnaClients
  def environment
    Polyphemus.instance.environment
  end

  def token
    Polyphemus.instance.config(:polyphemus, environment)[:token]
  end

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new(token: token, host: Polyphemus.instance.config(:magma, environment)[:host])
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new(token: token, host: Polyphemus.instance.config(:metis, environment)[:host])
  end
end

module WithLogger
  def logger
    Polyphemus.instance.logger
  end
end

WithEtnaClientsByEnvironment = EnvironmentScoped.new do
  include WithEtnaClients
end
