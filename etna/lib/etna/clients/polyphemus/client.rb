require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'

module Etna
  module Clients
    class Polyphemus
      def initialize(host:, token:, ignore_ssl: false)
        raise 'Polyphemus client configuration is missing host.' unless host
        raise 'Polyphemus client configuration is missing token.' unless token
        @etna_client = ::Etna::Client.new(
          host,
          token,
          routes_available: false,
          ignore_ssl: ignore_ssl)
      end

      def configuration(configuration_request = ConfigurationRequest.new)
        json = nil
        @etna_client.get(
          "/configuration",
          configuration_request) do |res|
          json = JSON.parse(res.body, symbolize_names: true)
        end

        ConfigurationResponse.new(json)
      end

      def job(job_request = JobRequest.new, &block)
        # Because this is a streaming response, just yield the response back.
        #   The consumer will have to iterate over the response.read_body, like
        #
        # polyphemus_client.job(request) do |response|
        #   response.read_body do |fragment|
        #     <fragment contains a chunk of text streamed back from the server>
        #   end
        # end
        @etna_client.post(
          "/#{job_request.project_name}/job",
          job_request) do |res|
            yield res
        end
      end
    end
  end
end
