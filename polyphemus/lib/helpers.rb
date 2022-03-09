require 'etna/clients'
require 'etna/environment_scoped'
require 'shellwords'

module WithSlackNotifications
  def notify_slack(message, messenger:self.class.name, channel:)
    `/bin/post-to-slack.sh #{Shellwords.escape(messenger)} #{Shellwords.escape(channel)} #{Shellwords.escape(message)} || true`
  end
end

module WithEtnaClients
  def environment
    Polyphemus.instance.environment
  end

  def token
    Polyphemus.instance.config(:polyphemus, environment)[:token]
  end

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new(
      token: token,
      host: Polyphemus.instance.config(:magma, environment)[:host])
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new(
      token: token,
      host: Polyphemus.instance.config(:metis, environment)[:host])
  end

  def metis_client=(client)
    @metis_client = client if Polyphemus.instance.test?
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
