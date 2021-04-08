require 'net/http/post/multipart'
require 'singleton'
require 'rack/utils'

module Etna
  class Client
    def initialize(host, token, routes_available: true, ignore_ssl: false)
      @host = host.sub(%r!/$!, '')
      @token = token
      @ignore_ssl = ignore_ssl

      if routes_available
        set_routes
        define_route_helpers
      end
    end

    attr_reader :routes

    def with_headers(headers, &block)
      @request_headers = headers.compact
      result = instance_eval(&block)
      @request_headers = nil
      return result
    end

    def signed_route_path(route, params)
      path = route_path(route,params)

      signatory = params.delete(:signatory)

      return path unless signatory

      hmac = Etna::Hmac.new(
        signatory,
        method: route[:method],
        host: URI(@host).host,
        path: path,
        expiration: (DateTime.now + 10).iso8601,
        id: signatory.id,
        nonce: SecureRandom.hex,
        headers: params.except(*route[:params].map(&:to_sym))
      )

      url_params = hmac.url_params(route[:method] == 'GET')

      return url_params[:path] + '?' + url_params[:query]
    end

    def route_path(route, params)
      Etna::Route.path(route[:route], params)
    end

    def multipart_post(endpoint, content, &block)
      uri = request_uri(endpoint)
      multipart = Net::HTTP::Post::Multipart.new uri.request_uri, content
      multipart.add_field('Authorization', "Etna #{@token}")
      request(uri, multipart, &block)
    end

    def post(endpoint, params = {}, &block)
      body_request(Net::HTTP::Post, endpoint, params, &block)
    end

    def get(endpoint, params = {}, &block)
      query_request(Net::HTTP::Get, endpoint, params, &block)
    end

    def head(endpoint, params = {}, &block)
      query_request(Net::HTTP::Head, endpoint, params, &block)
    end

    def options(endpoint, params = {}, &block)
      query_request(Net::HTTP::Options, endpoint, params, &block)
    end

    def delete(endpoint, params = {}, &block)
      body_request(Net::HTTP::Delete, endpoint, params, &block)
    end

    private

    def set_routes
      response = options('/')
      status_check!(response)
      @routes = JSON.parse(response.body, symbolize_names: true)
    end

    def define_route_helpers
      @routes.each do |route|
        next unless route[:name]
        self.define_singleton_method(route[:name]) do |params = {}|
          missing_params = (route[:params] - params.keys.map(&:to_s))

          unless missing_params.empty?
            raise ArgumentError, "Missing required #{missing_params.size > 1 ?
                'params' : 'param'} #{missing_params.join(', ')}"
          end

          response = send(
            route[:method].downcase,
            signed_route_path(route, params),
            params
          )

          if block_given?
            yield response
          else
            if response.content_type == 'application/json'
              return JSON.parse(response.body, symbolize_names: true)
            else
              return response.body
            end
          end
        end
      end
    end

    def body_request(type, endpoint, params = {}, &block)
      uri = request_uri(endpoint)
      req = type.new(uri.request_uri, request_headers)
      req.body = params.to_json
      request(uri, req, &block)
    end

    def query_request(type, endpoint, params = {}, &block)
      uri = request_uri(endpoint)

      if uri.query
        uri.query += "&" + URI.encode_www_form(params)
      else
        uri.query = URI.encode_www_form(params)
      end
      req = type.new(uri.request_uri, request_headers)
      request(uri, req, &block)
    end

    def request_uri(endpoint)
      URI("#{@host}#{endpoint}")
    end

    def request_headers
      {
          'Content-Type' => 'application/json',
          'Accept' => 'application/json, text/*',
          'Authorization' => "Etna #{@token}"
      }.update(
        @request_headers || {}
      )
    end

    def status_check!(response)
      status = response.code.to_i
      if status >= 400
        msg = response.content_type == 'application/json' ?
            json_error(response.body) :
            response.body
        raise Etna::Error.new(msg, status)
      end
    end

    def json_error(body)
      msg = JSON.parse(body, symbolize_names: true)
      if (msg.has_key?(:errors) && msg[:errors].is_a?(Array))
        return JSON.generate(msg[:errors])
      elsif msg.has_key?(:error)
        return JSON.generate(msg[:error])
      end
    end

    def request(uri, data)
      if block_given?
        verify_mode = @ignore_ssl ?
          OpenSSL::SSL::VERIFY_NONE :
          OpenSSL::SSL::VERIFY_PEER
        Net::HTTP.start(uri.host, uri.port, use_ssl: true, verify_mode: verify_mode) do |http|
          http.request(data) do |response|
            status_check!(response)
            yield response
          end
        end
      else
        verify_mode = @ignore_ssl ?
          OpenSSL::SSL::VERIFY_NONE :
          OpenSSL::SSL::VERIFY_PEER
        Net::HTTP.start(uri.host, uri.port, use_ssl: true, verify_mode: verify_mode) do |http|
          response = http.request(data)
          status_check!(response)
          return response
        end
      end
    end
  end
end
