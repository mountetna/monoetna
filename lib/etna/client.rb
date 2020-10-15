require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require 'rack/utils'

module Etna
  class Client
    def initialize(host, token, routes_available: true, persistent: true, ignore_ssl: false)
      @host = host.sub(%r!/$!, '')
      @token = token
      @persistent = persistent
      @ignore_ssl = ignore_ssl

      if routes_available
        set_routes
        define_route_helpers
      end
    end

    attr_reader :routes

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

          response = send(route[:method].downcase, route_path(route, params), params)
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

    def persistent_connection
      @http ||= begin
        http = Net::HTTP::Persistent.new
        http.read_timeout = 3600
        http.verify_mode = OpenSSL::SSL::VERIFY_NONE if @ignore_ssl
        http
      end
    end


    def body_request(type, endpoint, params = {}, &block)
      uri = request_uri(endpoint)
      req = type.new(uri.request_uri, request_params)
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
      req = type.new(uri.request_uri, request_params)
      request(uri, req, &block)
    end

    def request_uri(endpoint)
      URI("#{@host}#{endpoint}")
    end

    def request_params
      {
          'Content-Type' => 'application/json',
          'Accept' => 'application/json, text/*',
          'Authorization' => "Etna #{@token}"
      }
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
        if @persistent
          persistent_connection.request(uri, data) do |response|
            status_check!(response)
            yield response
          end
        else
          Net::HTTP.start(uri.host, uri.port, use_ssl: true) do |http|
            http.request(data) do |response|
              status_check!(response)
              yield response
            end
          end
        end
      else
        if @persistent
          response = persistent_connection.request(uri, data)
          status_check!(response)
          return response
        else
          Net::HTTP.start(uri.host, uri.port, use_ssl: true) do |http|
            response = http.request(data)
            status_check!(response)
            return response
          end
        end
      end
    end
  end
end
