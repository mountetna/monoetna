require 'etna/clients'
require 'etna/environment_scoped'
require 'shellwords'

module WithEtnaClients
  def environment
    Vesta.instance.environment
  end

  def token
    Vesta.instance.config(:vesta, environment)[:token]
  end

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new(
      token: token,
      host: Vesta.instance.config(:magma, environment)[:host])
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new(
      token: token,
      host: Vesta.instance.config(:metis, environment)[:host])
  end

  def metis_client=(client)
    @metis_client = client if Vesta.instance.test?
  end

  def janus_client
    @magma_client ||= Etna::Clients::Janus.new(
      token: token,
      host: Vesta.instance.config(:janus, environment)[:host])
  end
end

module WithLogger
  def logger
    Vesta.instance.logger
  end
end

WithEtnaClientsByEnvironment = EnvironmentScoped.new do
  include WithEtnaClients
end
