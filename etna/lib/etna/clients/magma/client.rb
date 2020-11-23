require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Magma
      attr_reader :host, :token, :ignore_ssl
      def initialize(host:, token:, persistent: true, ignore_ssl: false)
        raise 'Magma client configuration is missing host.' unless host
        raise 'Magma client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(
          host,
          token,
          routes_available: false,
          persistent: persistent,
          ignore_ssl: ignore_ssl)
        @host = host
        @token = token
        @ignore_ssl = ignore_ssl
      end

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
