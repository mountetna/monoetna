# a module to use for base application classes
# it is an included module rather than a class to spare
# us the burden of writing class Blah < Etna::Application
# whenever we want to use it as a container

require_relative './sign_service'
require_relative './command'
require_relative './generate_autocompletion_script'
require 'singleton'
require 'rollbar'
require 'yabeda'

module Etna::Application
  def self.included(other)
    other.include Singleton
    other.include Etna::CommandExecutor
    @@application = other

    other.const_set(:GenerateCompletionScript, Class.new(Etna::GenerateCompletionScript))
  end

  def self.find(klass)
    namespace = klass.name.split('::').first
    if (namespace_klass = Kernel.const_get(namespace)) && (namespace_klass.respond_to? :instance)
      return namespace_klass.instance
    end

    if @@application
      return @@application.instance
    end

    raise "Could not find application instance from #{namespace}, and not subclass of Application found."
  end

  # Used to find the application in development recorded vcr tests.
  # see spec/vcr.rb
  def dev_route
    "#{self.class.name.split('::').first.downcase}.development.local"
  end

  def self.register(app)
    @instance = app
  end

  def self.instance
    @instance
  end

  def initialize
    Etna::Application.register(self)
  end

  def configure(opts)
    @config = opts

    if (rollbar_config = config(:rollbar)) && rollbar_config[:access_token]
      Rollbar.configure do |config|
        config.access_token = rollbar_config[:access_token]
      end
    end
  end

  def setup_yabeda
    Yabeda.configure!
  end

  def write_job_metrics(name)
    node_metrics_dir = config(:node_metrics_dir) || "/tmp/metrics.prom"
    tmp_file = ::File.join(node_metrics_dir, "#{name}.prom.$$")
    ::File.open(tmp_file, "w") do |f|
      f.write(Prometheus::Client::Formats::Text.marshal(Prometheus::Client.registry))
    end

    require 'fileutils'
    ::FileUtils.mv(tmp_file, ::File.join(node_metrics_dir, "#{name}.prom"))
  end

  def setup_logger
    @logger = Etna::Logger.new(
      # The name of the log_file, required.
      config(:log_file),
      # Number of old copies of the log to keep.
      config(:log_copies) || 5,
      # How large the log can get before overturning.
      config(:log_size) || 1048576,
    )
    log_level = (config(:log_level) || 'warn').upcase.to_sym
    @logger.level = Logger.const_defined?(log_level) ? Logger.const_get(log_level) : Logger::WARN
  end

  # the application logger is available globally
  attr_reader :logger

  def config(type, env = environment)
    return nil if @config.nil?
    return nil if @config[env].nil?
    return nil unless @config[env].is_a?(Hash)
    @config[env][type]
  end

  def env_config(env = environment)
    return nil if @config.nil?
    return nil if @config[env].nil?
    return nil unless @config[env].is_a?(Hash)
    @config[env]
  end

  def sign
    @sign ||= Etna::SignService.new(self)
  end

  def environment
    (ENV["#{self.class.name.upcase}_ENV"] || :development).to_sym
  end

  def id
    ENV["APP_NAME"] || self.class.name.snake_case.split(/::/).last
  end

  def find_descendents(klass)
    ObjectSpace.each_object(Class).select do |k|
      k < klass
    end
  end

  def run_command(config, *args, &block)
    cmd, cmd_args, cmd_kwds = find_command(*args)
    cmd.setup(config)

    if block_given?
      return unless yield [cmd, cmd_args]
    end

    cmd.execute(*cmd.fill_in_missing_params(cmd_args), **cmd_kwds)
  rescue => e
    Rollbar.error(e)
    raise
  end
end

