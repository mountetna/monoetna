require_relative './etna/clients'

class EnvironmentScoped < Module
  def initialize(&block)
    environment_class = Class.new do
      class_eval(&block)

      attr_reader :environment
      def initialize(environment)
        @environment = environment
      end
    end

    super() do
      define_method :environment do |env|
        env = env.to_sym
        (@envs ||= {})[env] ||= environment_class.new(env)
      end
    end
  end
end

module WithEtnaClients
  def environment
    EtnaApp.instance.environment
  end

  def token
    if !ENV['TOKEN']
      puts "No environment variable TOKEN is set. Set your token with `export TOKEN=<your.janus.token>` before running Etna commands."
      exit
    end
    ENV['TOKEN']
  end

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new(token: token, host: EtnaApp.instance.config(:magma, environment)[:host])
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new(token: token, host: EtnaApp.instance.config(:metis, environment)[:host])
  end

  def janus_client
    @janus_client ||= Etna::Clients::Janus.new(token: token, host: EtnaApp.instance.config(:janus, environment)[:host])
  end
end

module WithLogger
  def logger
    EtnaApp.instance.logger
  end
end

WithEtnaClientsByEnvironment = EnvironmentScoped.new do
  include WithEtnaClients
end
