require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'
require_relative '../base_client'

module Etna
  module Clients
    class Janus < Etna::Clients::BaseClient
      def get_project(get_project_request = GetProjectRequest.new)
        json = nil
        @etna_client.get(
          "/api/admin/#{get_project_request.project_name}/info") do |res|
            json = JSON.parse(res.body, symbolize_names: true)
        end

        GetProjectResponse.new(json)
      end

      def get_projects()
        json = nil
        @etna_client.get('/api/user/projects') do |res|
          json = JSON.parse(res.body, symbolize_names: true)
        end

        GetProjectsResponse.new(json)
      end

      def add_project(add_project_request = AddProjectRequest.new)
        @etna_client.post('/api/admin/add_project', add_project_request) do |res|
          # Redirect, no response data
        end
      end

      def add_user(add_user_request = AddUserRequest.new)
        @etna_client.post(
          "/api/admin/#{add_user_request.project_name}/add_user",
          add_user_request) do |res|
          # Redirect, no response data
        end
      end

      def update_permission(update_permission_request = UpdatePermissionRequest.new)
        @etna_client.post(
          "/api/admin/#{update_permission_request.project_name}/update_permission",
          update_permission_request) do |res|
          # Redirect, no response data
        end
      end

      def refresh_token(refresh_token_request = RefreshTokenRequest.new)
        token = nil
        @etna_client.post('/api/tokens/generate') do |res|
          token = res.body
        end

        TokenResponse.new(token)
      end

      def validate_task_token
        @etna_client.post('/api/tokens/validate_task')
      end

      def get_nonce
        @etna_client.get('/api/tokens/nonce').body
      end

      def generate_token(token_type, signed_nonce: nil, project_name: nil, read_only: false)
        response = @etna_client.with_headers(
          'Authorization' => signed_nonce ? "Signed-Nonce #{signed_nonce}" : nil
        ) do
          post('/api/tokens/generate', token_type: token_type, project_name: project_name, read_only: read_only)
        end

        response.body
      end

      def get_project_stats(get_stats_request = GetStatsRequest.new)
        query = ""
        unless get_stats_request.project_names.nil?
          query = "?"
          query += get_stats_request.project_names.map do |name|
            "projects[]=#{name}"
          end.join('&')
        end

        json = nil

        @etna_client.get("/api/stats/projects#{query}") do |res|
          json = JSON.parse(res.body, symbolize_names: true)
        end

        json
      end
    end
  end
end
