require 'logger'

module Etna
  class Error < StandardError
    attr_reader :level, :status
    def initialize(msg = 'The application had an error')
      @level = Logger::WARN
      super
    end
  end

  class BadRequest < Etna::Error
    def initialize(msg = 'Client error', status = 422)
      super(msg)
      @status = status
    end
  end

  class ServerError < Etna::Error
    def initialize(msg = 'Server error', status = 500)
      super(msg)
      @status = status
      @level = Logger::ERROR
    end
  end
end
