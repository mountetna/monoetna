require 'net/http/post/multipart'
require 'singleton'
require_relative '../base_client'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class LocalMagmaClient # May only be used from within magma app.
      def logger
        ::Magma.instance.logger
      end

      def user
        permissions = Etna::Permissions.new([])
        permissions.add_permission(Etna::Permission.new('A', 'administration'))
        @user ||= Etna::User.new({name: 'Admin', email: 'etnaagent@ucsf.edu', perm: permissions.to_string})
      end

      def fake_request(request)
        Rack::Request.new({
          'rack.request.params' => request.as_json,
          'etna.logger' => self.logger,
          'etna.user' => user,
          'etna.request_id' => 'local',
        })
      end

      def parse_json_response(response)
        status, headers, body = response
        body = body.join('')
        if status < 200 || status >= 300
          raise Etna::Error.new(body, status)
        end

        JSON.parse(body)
      end

      def retrieve(retrieval_request = RetrievalRequest.new)
        controller = ::RetrieveController.new(fake_request(retrieval_request), nil)
        Magma::RetrievalResponse.new(parse_json_response(controller.action))
      end

      def update_model(update_model_request = UpdateModelRequest.new)
        controller = ::UpdateModelController.new(fake_request(update_model_request), nil)
        Magma::UpdateModelResponse.new(parse_json_response(controller.action))
      end
    end

    class Magma < Etna::Clients::BaseClient
      # This endpoint returns models and records by name:
      # e.g. params:
      # {
      #   model_name: "model_one", # or "all"
      #   record_names: [ "rn1", "rn2" ], # or "all",
      #   attribute_names:  "all"
      # }
      def retrieve(retrieval_request = RetrievalRequest.new)
        json = nil
        @etna_client.post('/retrieve', retrieval_request) do |res|
          json = JSON.parse(res.body)
        end

        RetrievalResponse.new(json)
      end

      # This 'query' end point is used to fetch data by graph query
      # See question.rb for more detail
      def query(query_request = QueryRequest.new)
        json = nil
        @etna_client.post('/query', query_request) do |res|
          json = JSON.parse(res.body)
        end

        QueryResponse.new(json)
      end

      def update(update_request = UpdateRequest.new)
        json = nil
        @etna_client.multipart_post('/update', update_request.encode_multipart_content) do |res|
          json = JSON.parse(res.body)
        end

        UpdateResponse.new(json)
      end

      def update_json(update_request = UpdateRequest.new)
        json = nil
        @etna_client.post('/update', update_request) do |res|
          json = JSON.parse(res.body)
        end

        UpdateResponse.new(json)
      end

      def update_model(update_model_request = UpdateModelRequest.new)
        json = nil
        @etna_client.post('/update_model', update_model_request) do |res|
          json = JSON.parse(res.body)
        end

        UpdateModelResponse.new(json)
      end
    end
  end
end
