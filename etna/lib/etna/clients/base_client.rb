require 'base64'
require 'json'
require 'date'

module Etna
  module Clients
    class BaseClient
      attr_reader :host, :token, :ignore_ssl
      def initialize(host:, token:, ignore_ssl: false)
        raise "#{self.class.name} client configuration is missing host." unless host
        raise "#{self.class.name} client configuration is missing token." unless token

        @token = token
        raise "Your token is expired." if token_expired?

        @etna_client = ::Etna::Client.new(
          host,
          token,
          routes_available: false,
          ignore_ssl: ignore_ssl)
        @host = host
        @ignore_ssl = ignore_ssl
      end

      def token_expired?
        # Has the token already expired?
        token_will_expire?(0)
      end

      def token_will_expire?(offset=3000)
        # offset in seconds
        # Will the user's token expire in the given amount of time?
        epoch_seconds = JSON.parse(Base64.urlsafe_decode64(token.split('.')[1]))["exp"]
        expiration = DateTime.strptime(epoch_seconds.to_s, "%s")
        expiration <= DateTime.now.new_offset + offset
      end
    end
  end
end
