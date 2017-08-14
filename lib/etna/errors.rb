require 'logger'

module Etna
  class Error < StandardError
    attr_reader :level
    def initialize(msg="The application had an error")
      @level = Logger::WARN
      super
    end
  end
end
