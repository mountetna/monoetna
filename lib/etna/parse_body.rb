module Etna
  class ParseBody
    def initialize(app)
      @app = app
    end

    def call(env)
      params = env['rack.request.params'] || {}

      case env['CONTENT_TYPE']
      when %r{application/json}i
        body = env['rack.input'].read
        if body =~ %r/^\s*\{/
          params.update(
            JSON.parse(body)
          )
        end
      when %r{application/x-www-form-urlencoded}i
        params.update(
          Rack::Utils.parse_nested_query(
            env['rack.input'].read
          )
        )
      when %r{multipart/form-data}i
        params.update(
          Rack::Multipart.parse_multipart(env)
        )
      end
      # Always parse the params that are url-encoded.
      params.update(
        Rack::Utils.parse_nested_query(
          env['QUERY_STRING'], '&'
        )
      )
      env.update(
        'rack.request.params' => params
      )
      @app.call(env)
    end
  end
end
