require 'net/http/post/multipart'
require 'singleton'
require 'rack/utils'

module Etna
  module ClientUtils
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
  end
  class Client
    include Etna::ClientUtils



    def initialize(host, token, routes_available: true, ignore_ssl: false, max_retries: 10, backoff_time: 15)
      @host = host.sub(%r!/$!, '')
      @token = token
      @ignore_ssl = ignore_ssl
      @max_retries = max_retries
      @backoff_time = backoff_time

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
      retrier.retry_request(uri, multipart, &block)
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

    def retrier
      @retrier ||= Retrier.new(ignore_ssl: @ignore_ssl, max_retries: @max_retries, backoff_time: @backoff_time)
    end

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
      retrier.retry_request(uri, req, &block)
    end

    def query_request(type, endpoint, params = {}, &block)
      uri = request_uri(endpoint)

      if uri.query
        uri.query += "&" + URI.encode_www_form(params)
      else
        uri.query = URI.encode_www_form(params)
      end
      req = type.new(uri.request_uri, request_headers)
      retrier.retry_request(uri, req, &block)
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

    class Retrier
      # Ideally the retry code would be centralized with metis_client ...
      #   unsure what would be the best approach to do that, at this moment.
      include Etna::ClientUtils

      def initialize(ignore_ssl: false, max_retries: 10, backoff_time: 15)
        @max_retries = max_retries
        @ignore_ssl = ignore_ssl
        @backoff_time = backoff_time
      end

      def retry_request(uri, data, retries: 0, &block)
        retries += 1

        begin
          logger.write("\rWaiting for server restart"+"."*retries+"\x1b[0K")

          sleep @backoff_time
        end if retries > 1

        if retries < @max_retries
          begin
            if block_given?
              request(uri, data, &block)
            else
              response = request(uri, data)
            end
          rescue OpenSSL::SSL::SSLError => e
            if e.message =~ /write client hello/
              logger.write("SSL error, retrying")
              return retry_request(uri, data, retries: retries)
            end
            raise e
          rescue *net_exceptions => e
            logger.write("Received #{e.class.name}, retrying")
            return retry_request(uri, data, retries: retries)
          end

          begin
            retry_codes = ['503', '502', '504', '408']
            if retry_codes.include?(response.code)
              logger.write("Received response with code #{response.code}, retrying")
              return retry_request(uri, data, retries: retries)
            elsif response.code == '500' && response.body.start_with?("Puma caught")
              logger.write("Received 500 Puma error #{response.body.split("\n").first}, retrying")
              return retry_request(uri, data, retries: retries)
            end

            return response
          end unless block_given?
        end

        raise ::Etna::Error, "Could not contact server, giving up" unless block_given?
      end

      private

      def request(uri, data)
        if block_given?
          verify_mode = @ignore_ssl ?
            OpenSSL::SSL::VERIFY_NONE :
            OpenSSL::SSL::VERIFY_PEER
          Net::HTTP.start(uri.host, uri.port, use_ssl: true, verify_mode: verify_mode, read_timeout: 300) do |http|
            http.request(data) do |response|
              status_check!(response)
              yield response
            end
          end
        else
          verify_mode = @ignore_ssl ?
            OpenSSL::SSL::VERIFY_NONE :
            OpenSSL::SSL::VERIFY_PEER
          Net::HTTP.start(uri.host, uri.port, use_ssl: true, verify_mode: verify_mode, read_timeout: 300) do |http|
            response = http.request(data)
            status_check!(response)
            return response
          end
        end
      end

      def net_exceptions
        retry_exceptions = [
          Errno::ECONNREFUSED,
          Errno::ECONNRESET,
          Errno::ENETRESET,
          Errno::EPIPE,
          Errno::ECONNABORTED,
          Errno::EHOSTDOWN,
          Errno::EHOSTUNREACH,
          Errno::EINVAL,
          Errno::ETIMEDOUT,
          Net::ReadTimeout,
          Net::HTTPFatalError,
          Net::HTTPBadResponse,
          Net::HTTPHeaderSyntaxError,
          Net::ProtocolError,
          Net::HTTPRequestTimeOut,
          Net::HTTPGatewayTimeOut,
          Net::HTTPBadRequest,
          Net::HTTPBadGateway,
          Net::HTTPError,
          Net::HTTPInternalServerError,
          Net::HTTPRetriableError,
          Net::HTTPServerError,
          Net::HTTPServiceUnavailable,
          Net::HTTPUnprocessableEntity,
          Net::OpenTimeout,
          IOError,
          EOFError,
          Timeout::Error
        ]

        begin
          retry_exceptions << Net::HTTPRequestTimeout
          retry_exceptions << Net::HTTPGatewayTimeout
          retry_exceptions << Net::WriteTimeout
        end if RUBY_VERSION > "2.5.8"

        retry_exceptions
      end
    end
  end
end
