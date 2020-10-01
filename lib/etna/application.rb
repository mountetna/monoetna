# a module to use for base application classes
# it is an included module rather than a class to spare
# us the burden of writing class Blah < Etna::Application
# whenever we want to use it as a container

require_relative './sign_service'
require 'singleton'
require 'rollbar'

module Etna::Application
  def self.included(other)
    other.include Singleton
  end

  def self.find(klass)
    Kernel.const_get(
      klass.name.split('::').first
    ).instance
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

  # In some cases, such as utility scripts that span across environments, it may be necessary to override
  # the environment source.
  def config(type, env = environment)
    @config[env][type]
  end

  def sign
    @sign ||= Etna::SignService.new(self)
  end

  def environment
    (ENV["#{self.class.name.upcase}_ENV"] || :development).to_sym
  end

  def find_descendents(klass)
    ObjectSpace.each_object(Class).select do |k|
      k < klass
    end
  end

  def run_command(config, cmd = :help, *args)
    cmd = cmd.to_sym
    if commands.key?(cmd)
      commands[cmd].setup(config)
      commands[cmd].execute(*args)
    else
      commands[:help].execute
    end
  rescue => e
    Rollbar.error(e)
    raise
  end

  def commands
    @commands ||= Hash[
      find_descendents(Etna::Command).map do |c|
        cmd = c.new
        [ cmd.name, cmd ]
      end
    ]
  end
end

