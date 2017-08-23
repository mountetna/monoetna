require 'logger'

module Etna
  class Error < StandardError
    attr_reader :level
    def initialize(msg='The application had an error')
      @level = Logger::WARN
      super
    end
  end

  class BadRequest < Etna::Error
  end

  class ServerError < Etna::Error
    def initialize(msg='Server error')
      super
      @level = Logger::ERROR
    end
  end
end
