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

      def get_config(project_name, config_id, version_number)
        json = nil
        @etna_client.post("/api/etl/#{project_name}/configs/#{config_id}", version_number: version_number) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def get_runtime_config(project_name, config_id)
        json = nil
        @etna_client.get("/api/etl/#{project_name}/runtime_configs/#{config_id}") do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def update_run(project_name, run_id, updates)
        payload = {
            run_id: run_id,
            config_id: updates[:config_id],
            version_number: updates[:version_number],
            state: updates[:state],
            orchestrator_metadata: updates[:orchestrator_metadata],
            output: updates[:output],
            append_output: updates[:append_output]
        }.compact          

        json = nil
        @etna_client.post("/api/etl/#{project_name}/run/update/#{run_id}", payload) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def get_run(project_name, run_id)
        json = nil
        @etna_client.get("/api/etl/#{project_name}/run/#{run_id}") do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def get_previous_run(project_name, config_id, version_number)
        json = nil
        puts "Requesting in client: /api/etl/#{project_name}/run/previous/#{config_id}"
        @etna_client.post("/api/etl/#{project_name}/run/previous/#{config_id}", version_number: version_number) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def update_runtime_config(project_name, config_id, updates)
        payload = {
          config_id: updates[:config_id],
          runtime_config: updates[:runtime_config],
          run_interval: updates[:run_interval]
        }.compact          

        json = nil
        @etna_client.post("/api/etl/#{project_name}/runtime_config/#{config_id}", payload) do |res|
          json = JSON.parse(res.body)
        end
        json
      end

      def get_run_metadata(project_name, config_id)
        json = nil
        @etna_client.get("/api/etl/#{project_name}/runtime_config/#{config_id}") do |res|
          json = JSON.parse(res.body)
        end
        json
      end

    end
  end
end
