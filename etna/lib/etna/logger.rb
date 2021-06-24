require 'logger'
require 'rollbar'

module Etna
  class Logger < ::Logger
    def initialize(log_dev, age, size=1048576)
      # On windows, these posix devices exist, but are not mounted in *nix style paths.
      # Swap the paths out with the actual IO handles instead.
      if log_dev == '/dev/stdout'
        log_dev = STDOUT
      elsif log_dev == '/dev/stderr'
        log_dev = STDERR
      end

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

      if defined? Yabeda
        Yabeda.etna.rollbar_errors.increment({}, 1) rescue nil
      end

      Rollbar.error(e)
    end

    def log_request(request)
      request.env['etna.logger'] = self
      request.env['etna.request_id'] = (rand*36**6).to_i.to_s(36)
    end
  end
end
