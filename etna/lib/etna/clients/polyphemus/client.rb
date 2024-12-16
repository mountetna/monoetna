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

      def get_config(config_id, version_number)
        json = nil
        @etna_client.post("/api/workflows/#{@project_name}/configs/#{config_id}", version_number: version_number) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def update_workflow_run(run_id, state)
        json = nil
        @etna_client.post("/api/workflows/#{@project_name}/workflow_state/update/#{config_id}", {
          run_id: run_id,
          state: state
        }) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def get_workflow_run(run_id)
        json = nil
        @etna_client.get("/api/workflows/#{@project_name}/run/#{run_id}") do |res|
          json = JSON.parse(res.body)
        end
        json[:state]
      end

      def update_run_metadata(run_id, workflow_json)
        def update_run_metadata(run_id, workflow_json)
          payload = {
            run_id: run_id,
            config_id: workflow_json[:config_id],
            version_number: workflow_json[:version_number],
            orchestrator_metadata: workflow_json[:orchestrator_metadata],
            runtime_config: workflow_json[:runtime_config],
            output: workflow_json[:output],
            append_output: workflow_json[:append_output],
            run_interval: workflow_json[:run_interval]
          }.compact  # Removes any key-value pairs where the value is nil
        
          # Make the POST request with the cleaned payload
          json = nil
          @etna_client.post("/api/workflows/#{@project_name}/run_metadata", payload) do |res|
            json = JSON.parse(res.body)
          end
          json
        end


      end

      def get_run_metadata(run_id)
        json = nil
        @etna_client.get("/api/workflows/#{@project_name}/run_metadata/#{run_id}") do |res|
          json = JSON.parse(res.body)
        end
        json
      end

    end
  end
end
