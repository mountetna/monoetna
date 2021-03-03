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
