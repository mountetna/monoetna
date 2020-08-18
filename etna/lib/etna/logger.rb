require 'logger'
require 'rollbar'

module Etna
  class Logger < ::Logger
    def initialize(log_dev, age, size, rollbar:)
      super

      self.formatter = proc do |severity, datetime, progname, msg|
        format(severity, datetime, progname, msg)
      end

      if rollbar && rollbar[:access_token]
        Rollbar.configure do |config|
          config.access_token = rollbar[:access_token]
        end
        @configured_rollbar = true
      end
    end

    def format(severity, datetime, progname, msg)
      "#{severity}:#{datetime.iso8601} #{msg}\n"
    end

    def log_error(e)
      error(e.message)
      e.backtrace.each do |trace|
        error(trace)
      end

      if @configured_rollbar
        Rollbar.error(e)
      end
    end

    def log_request(request)
      request.env['etna.logger'] = self
      request.env['etna.request_id'] = (rand*36**6).to_i.to_s(36)
    end
  end
end
