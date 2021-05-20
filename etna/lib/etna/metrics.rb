
module Etna
  class MetricsExporter
    def initialize(app, path: '/metrics')
      @app = app
      @path = path
    end

    def exporter
      @exporter ||= begin
        exporter = Yabeda::Prometheus::Exporter.new(@app, path: @path)
        Rack::Auth::Basic.new(exporter) do |user, pw|
          user == 'prometheus' && pw == ENV['METRICS_PW']
        end
      end
    end

    def call(env)
      if env['PATH_INFO'] == @path
        exporter.call(env)
      else
        @app.call(env)
      end
    end
  end
end