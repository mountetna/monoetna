# a module to use for base application classes
# it is an included module rather than a class to spare
# us the burden of writing class Blah < Etna::Application
# whenever we want to use it as a container

require_relative './sign_service'
require 'singleton'

module Etna::Application
  def self.included(other)
    other.include Singleton
  end

  def self.find(klass)
    Kernel.const_get(
      klass.name.split('::').first
    ).instance
  end

  def configure(opts)
    @config = opts
  end

  def config(type)
    @config[environment][type]
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

