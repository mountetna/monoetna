require 'net/http/persistent'
require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    module Magma
      class Client
        class Client
          def initialize(host, token)
            @etna_client = ::Etna::Client.new(host, token, routes_available: false)
            raise 'Magma client configuration is missing host.' unless @host
            raise 'Magma client configuration is missing token.' unless @token
          end

          # This endpoint returns models and records by name:
          # e.g. params:
          # {
          #   model_name: "model_one", # or "all"
          #   record_names: [ "rn1", "rn2" ], # or "all",
          #   attribute_names:  "all"
          # }
          def retrieve(retrieval_request)
            json = @etna_client.post('/retrieve', retrieval_request) do |res|
              JSON.parse(res.body)
            end

            RetrievalResponse.new(json)
          end

          # This 'query' end point is used to fetch data by graph query
          # See question.rb for more detail
          def query(query_request)
            json = @etna_client.post('/query', query_request) do |res|
              JSON.parse(res.body)
            end

            QueryResponse.new(json)
          end

          # Post revisions to Magma records
          # { model_name: { record_name: { attribute1: 1, attribute2: 2 } } } }
          # data can also be a File or IO stream
          def update(update_request)
            json = @etna_client.multipart_post('/update', update_request.encode_multipart_content) do |res|
              JSON.parse(res.body)
            end

            UpdateResponse.new(json)
          end

          def update_model(update_model_request)
            json = @etna_client.post('/update_model', update_model_request) do |res|
              JSON.parse(res.body)
            end

            UpdateModelResponse.new(json)
          end

          private

          def persistent_connection
            @http ||= begin
                        http = Net::HTTP::Persistent.new
                        http.read_timeout = 3600
                        http
                      end
          end

          def json_post(endpoint, token, params, status_errors = {}, &block)
            post(endpoint, 'application/json', token, params.to_json, status_errors, &block)
          end

          def multipart_post(endpoint, token, content, status_errors = {})
            uri = URI("#{@host}/#{endpoint}")
            multipart = Net::HTTP::Post::Multipart.new uri.path, content
            multipart.add_field('Authorization', "Etna #{token}")

            request(uri, multipart, status_errors)
          end

          def post(endpoint, content_type, token, body, status_errors, &block)
            uri = URI("#{@host}/#{endpoint}")
            post = Net::HTTP::Post.new(
                uri.path,
                'Content-Type' => content_type,
                'Accept' => 'application/json, text/*',
                'Authorization' => "Etna #{token}"
            )
            post.body = body
            request(uri, post, status_errors, &block)
          end

          def status_check(response, status_errors)
            status = response.code.to_i
            if status >= 500
              raise Magma::ClientError.new(status, (status_errors[500] || {}).merge(errors: ['A Magma server error occured.']))
            elsif status >= 400
              errors = response.content_type == 'application/json' ? JSON.parse(response.body)['errors'] : [response.body]
              raise Magma::ClientError.new(status, (status_errors[400] || {}).merge(errors: errors))
            end
          end

          def request(uri, data, status_errors, &block)
            if block_given?
              persistent_connection.request(uri, data) do |response|
                status_check(response, status_errors)
                yield response
              end
            else
              response = persistent_connection.request(uri, data)
              status_check(response, status_errors)
              return response
            end
          end
        end
      end
    end
  end
end
