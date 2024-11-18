require 'net/http/post/multipart'
require 'singleton'
require_relative '../../client'
require_relative './models'
require_relative '../base_client'

module Etna
  module Clients
    class Polyphemus < Etna::Clients::BaseClient
      def configuration(configuration_request = ConfigurationRequest.new)
        json = nil
        @etna_client.get(
          "/configuration",
          configuration_request) do |res|
          json = JSON.parse(res.body, symbolize_names: true)
        end

        ConfigurationResponse.new(json)
      end

      def job(job_request = JobRequest.new)
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

      def get_workflow(project_name, workflow_name, revision: "latest")
      end

      def get_workflow_state(argo_id)
      end

      def update_workflow_state(argo_id, state)
      end

    end
  end
end
