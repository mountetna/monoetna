module Etna
  class ParseBody
    def initialize(app)
      @app = app
    end

    def call(env)
      case env['CONTENT_TYPE']
      when %r{application/json}i
        body = env['rack.input'].read
        if body =~ %r/^\s*\{/
          env.update(
            'rack.request.params'=> JSON.parse(body)
          )
        end
      when %r{application/x-www-form-urlencoded}i
        env.update(
          'rack.request.params'=> Rack::Utils.parse_nested_query(
            env['rack.input'].read
          )
        )
      when %r{multipart/form-data}i
        env.update(
          'rack.request.params'=> Rack::Multipart.parse_multipart(env)
        )
      else
        # Assume the params are url-encoded.
        env.update(
          'rack.request.params'=> Rack::Utils.parse_nested_query(
            env['QUERY_STRING'], '&'
          )
        )
      end
      @app.call(env)
    end
  end
end
