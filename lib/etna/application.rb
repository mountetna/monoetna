# a module to use for base application classes
# it is an included module rather than a class to spare
# us the burden of writing class Blah < Etna::Application
# whenever we want to use it as a container

require 'singleton'

module Etna::Application
  def self.included(other)
    other.include Singleton
  end

  def configure(opts)
    @config = opts
  end

  def config(type)
    @config[environment][type]
  end

  def environment
    (ENV["#{self.class.name.upcase}_ENV"] || :development).to_sym
  end

  def find_descendents(klass)
    ObjectSpace.each_object(Class).select do |k|
      k < klass
    end
  end
end

