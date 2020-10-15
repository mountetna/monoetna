require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Polyphemus
      def initialize(host:, token:, ignore_ssl: false, persistent: true)
        raise 'Polyphemus client configuration is missing host.' unless host
        raise 'Polyphemus client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(
          host,
          token,
          routes_available: false,
          persistent: persistent,
          ignore_ssl: ignore_ssl)
      end

      def configuration(configuration_request = ConfigurationRequest.new)
        json = nil
        @etna_client.get(
          "/configuration",
          configuration_request) do |res|
          json = JSON.parse(res.body)
        end

        ConfigurationResponse.new(json)
      end
    end
  end
end
