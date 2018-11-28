require 'logger'

module Etna
  class Logger < ::Logger
    def initialize(log_dev, age, size)
      super

      self.formatter = proc do |severity, datetime, progname, msg|
        format(severity, datetime, progname, msg)
      end
    end

    def format(severity, datetime, progname, msg)
      "#{severity} #{datetime.iso8601} #{@request.env['etna.request_id']} #{msg}\n"
    end

    def request=(request)
      @request = request
      @request.env['etna.logger'] = self
      @request.env['etna.request_id'] = (rand*36**6).to_i.to_s(36)
    end
  end
end
