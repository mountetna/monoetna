require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'
require_relative '../base_client'

module Etna
  module Clients
    class Janus < Etna::Clients::BaseClient
      def get_project(get_project_request = GetProjectRequest.new)
        html = nil
        @etna_client.get(
          "/project/#{get_project_request.project_name}",
          get_project_request) do |res|
          html = res.body
        end

        HtmlResponse.new(html)
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

      def update_permission(update_permission_request = UpdatePermissionRequest.new)
        @etna_client.post(
          "/update_permission/#{update_permission_request.project_name}",
          update_permission_request) do |res|
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

      def viewer_token(viewer_token_request = ViewerTokenRequest.new)
        token = nil
        @etna_client.get('/viewer_token', viewer_token_request) do |res|
          token = res.body
        end

        TokenResponse.new(token)
      end

      def validate_task_token(validate_task_token_request = ValidateTaskTokenRequest.new)
        token = nil
        @etna_client.post('/api/tokens/task/validate', validate_task_token_request)
      end
    end
  end
end
