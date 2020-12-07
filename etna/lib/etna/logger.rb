require 'logger'
require 'rollbar'

module Etna
  class Logger < ::Logger
    def initialize(log_dev, age, size)
      super
      self.formatter = proc do |severity, datetime, progname, msg|
        format(severity, datetime, progname, msg)
      end
    end

    def initialize(log_dev, age)
      super
      self.formatter = proc do |severity, datetime, progname, msg|
        format(severity, datetime, progname, msg)
      end
    end

    def format(severity, datetime, progname, msg)
      "#{severity}:#{datetime.iso8601} #{msg}\n"
    end

    def warn(msg, &block)
      super
    end

    def error(msg, &block)
      super
    end

    def fatal(msg, &block)
      super
    end

    def log_error(e)
      error(e.message)
      e.backtrace.each do |trace|
        error(trace)
      end

      Rollbar.error(e)
    end

    def log_request(request)
      request.env['etna.logger'] = self
      request.env['etna.request_id'] = (rand*36**6).to_i.to_s(36)
    end
  end
end
