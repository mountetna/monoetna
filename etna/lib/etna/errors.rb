require 'logger'

module Etna
  class Error < StandardError
    attr_reader :level, :status
    def initialize(msg = 'The application had an error', status=500)
      super(msg)
      @level = Logger::WARN
      @status = status
    end
  end

  class Forbidden < Etna::Error
    def initialize(msg = 'Action not permitted', status = 403)
      super
    end
  end

  class Unauthorized < Etna::Error
    def initialize(msg = 'Unauthorized request', status = 401)
      super
    end
  end

  class NotFound < Etna::Error
    def initialize(msg = 'Resource not found', status = 404)
      super
    end
  end

  class BadRequest < Etna::Error
    def initialize(msg = 'Client error', status = 422)
      super
    end
  end

  class ServerError < Etna::Error
    def initialize(msg = 'Server error', status = 500)
      super
      @level = Logger::ERROR
    end
  end

  class TooManyRequests < Etna::Error
    def initialize(msg = 'Too many requests', status = 429)
      super
      @level = Logger::ERROR
    end
  end

end
