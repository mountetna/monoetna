require 'base64'
require 'json'
require 'date'

module Etna
  module Clients
    class BaseClient
      attr_reader :host, :token, :ignore_ssl
      def initialize(host:, token:, ignore_ssl: false, logger: nil, routes_available: false)
        raise "#{self.class.name} client configuration is missing host." unless host

        @token = token
        raise "Your token is expired." if token && token_expired?

        @etna_client = ::Etna::Client.new(
          host,
          token,
          routes_available: routes_available,
          ignore_ssl: ignore_ssl,
          logger: logger)
        @host = host
        @ignore_ssl = ignore_ssl
      end

      def token_expired?
        # Has the token already expired?
        token_will_expire?(0)
      end

      def safe_parse(response)
        if response['Content-Type'] == 'application/json'
          return JSON.parse(response.body)
        end

        raise "Could not parse non-JSON response #{response.code} #{response.body} #{response.each_header.to_h}"
      end

      def token_will_expire?(offset=3000)
        # offset in seconds
        # Will the user's token expire in the given amount of time?
        payload = JSON.parse(Base64.urlsafe_decode64(token.split('.')[1]))
        epoch_seconds = payload["exp"]
        expiration = DateTime.strptime(epoch_seconds.to_s, "%s")
        expiration <= DateTime.now.new_offset + offset
      end
    end
  end
end
