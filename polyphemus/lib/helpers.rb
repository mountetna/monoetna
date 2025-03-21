require 'etna/clients'
require 'etna/environment_scoped'
require 'shellwords'
require 'net/http'
require 'uri'

module WithSlackNotifications
  MESSAGE_SIZE=8000

  def chunk_message(message)
    return message.scan(/.{1,#{MESSAGE_SIZE}}/m)
  end

  def notify_slack(message, messenger: self.class.name, channel: nil, webhook_url: Polyphemus.instance.config(:slack_webhook_url))
    uri = URI(webhook_url)

    messages = chunk_message(message)

    Net::HTTP.start(uri.hostname, uri.port, use_ssl: true) do |http|
      messages.each do |msg|
        req = Net::HTTP::Post.new(uri, 'Content-Type' => 'application/json')

        req.body = {
          channel: channel,
          username: "Polyphemus",
          text: msg,
          icon_emoji: ":polyphemus:"
        }.compact.to_json

        http.request(req)
      end
    end
  end
end

module WithEtnaClients
  def environment
    Polyphemus.instance.environment
  end

  def token
    @token || Polyphemus.instance.config(:polyphemus, environment)[:token]
  end

  def janus_client
    @janus_client ||= Etna::Clients::Janus.new(
      token: token,
      host: Polyphemus.instance.config(:janus, environment)[:host])
  end

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new(
      token: token,
      host: Polyphemus.instance.config(:magma, environment)[:host])
  end

  def gnomon_client(logger: nil)
    @gnomon_client ||= Etna::Clients::Gnomon.new(
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

  def polyphemus_client
    @polyphemus_client ||= Etna::Clients::Polyphemus.new(
      token: token,
      host: Polyphemus.instance.config(:polyphemus, environment)[:host])
  end

  def reset_clients!
    @janus_client = nil
    @magma_client = nil
    @gnomon_client = nil
    @metis_client = nil
    @polyphemus_client = nil
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
