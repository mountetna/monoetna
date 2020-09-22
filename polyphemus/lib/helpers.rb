require 'etna/clients'

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
  def project
    raise "project must be implemented in subclasses!"
  end

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
