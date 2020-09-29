require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require_relative './models'

module Etna
  module Clients
    class Janus
      def initialize(host:, token:, persistent: true)
        raise 'Janus client configuration is missing host.' unless host
        raise 'Janus client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(host, token, routes_available: false, persistent: persistent)
      end

      def add_project(add_project_request = AddProjectRequest.new)
        @etna_client.post('/add_project', add_project_request) do |res|
          # Redirect, no response data
        end
      end

      def add_user(add_user_request = AddUserRequest.new)
        @etna_client.post(
          "/add_user/#{add_user_request.project_name}",
          add_user_request) do |res|
          # Redirect, no response data
        end
      end

      def refresh_token(refresh_token_request = RefreshTokenRequest.new)
        token = nil
        @etna_client.get('/refresh_token', refresh_token_request) do |res|
          token = res.body
        end

        TokenResponse.new(token)
      end
    end
  end
end
