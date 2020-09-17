require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require 'rack/utils'
require 'vcr'

module Etna
  class Client
    def initialize(host, token, routes_available: true)
      host = host.sub(%r!/$!,'')

      @host, @cassette_name = parse_cassette_from_host(host)

      if @cassette_name
        configure_vcr
      end

      @token = token

      if routes_available
        set_routes
        define_route_helpers
      end
    end

    attr_reader :routes

    def configure_vcr
      @@configured_vcr ||= begin
        VCR.configure do |c|
          c.hook_into :webmock
          c.cassette_library_dir = ::File.join(__dir__, '..', '..', 'spec', 'fixtures', 'cassettes')
          c.default_cassette_options = { record: ENV['IS_CI'] ? :none : :once }

          # Filter the authorization headers of any request by replacing any occurrence of that request's
          # Authorization value with <AUTHORIZATION>
          c.filter_sensitive_data('<AUTHORIZATION>') do |interaction|
            interaction.request.headers['Authorization']
          end
        end
      end
    end

    # Returns a canonical http URI for the host, and also returns the file in which a recording should be made / read from.
    def self.parse_cassette_from_host(host_string)
      uri = URI.parse(host_string)
      if uri.scheme == 'vcr'
        uri.scheme = 'https'
        cassette_name = '--unnamed--'
        if uri.query
          cassette_name = Rack::Utils.parse_nested_query(uri.query)['cassette'] || cassette_name
        end
        [uri.to_s, cassette_name]
      end

      [host_string, nil]
    end

    def route_path(route, params)
      Etna::Route.path(route[:route], params)
    end

    def multipart_post(endpoint, content, &block)
      uri = request_uri(endpoint)
      multipart = Net::HTTP::Post::Multipart.new uri.path, content
      multipart.add_field('Authorization', "Etna #{@token}")
      request(uri, multipart, &block)
    end

    def post(endpoint, params={}, &block)
      body_request(Net::HTTP::Post, endpoint, params, &block)
    end

    def get(endpoint, params={}, &block)
      query_request(Net::HTTP::Get, endpoint, params, &block)
    end

    def options(endpoint, params={}, &block)
      query_request(Net::HTTP::Options, endpoint, params, &block)
    end

    def delete(endpoint, params={}, &block)
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
        self.define_singleton_method(route[:name]) do |params={}|

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
                  http
                end
    end


    def body_request(type, endpoint, params={}, &block)
      uri = request_uri(endpoint)
      req = type.new(uri.request_uri,request_params)
      req.body = params.to_json
      request(uri, req, &block)
    end

    def query_request(type, endpoint, params={}, &block)
      uri = request_uri(endpoint)
      uri.query = URI.encode_www_form(params)
      req = type.new(uri.request_uri, request_params)
      request(uri, req, &block)
    end

    def request_uri(endpoint)
      URI("#{@host}#{endpoint}")
    end

    def request_params
      {
        'Content-Type' => 'application/json',
        'Accept'=> 'application/json, text/*',
        'Authorization'=>"Etna #{@token}"
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

    def request(uri, data, in_recording: false, &block)
      if !in_recording && @cassette_name
        VCR.use_cassette(@cassette_name) do
          return request(uri, data, in_recording: true, &block)
        end
      end

      if block_given?
        persistent_connection.request(uri, data) do |response|
          status_check!(response)
          yield response
        end
      else
        response = persistent_connection.request(uri, data)
        status_check!(response)
        return response
      end
    end
  end
end
