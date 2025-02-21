class Polyphemus
  class IngestionManifest < Polyphemus::WorkflowManifest
    class << self
      def as_json
        {
          name: 'ingestion',
          schema: {
            "$schema": "http://json-schema.org/draft-07/schema#",
            type: "object",
            description: "This loader downloads data from a remote host to Metis or another remote host",
	    title: "Ingestion Loader",
            properties: {
              bucket_name: { type: 'string' },
              magic_string: { type: 'string' },
              ingest_root_path: { type: 'string' },
              metis_root_path: { type: 'string' },
              deposit_root_path: { type: 'string' }
            }
          },
          runtime_params: {
            initial_start_scan_time: {
              type: 'integer', # unix timestamp
              description: 'Scan start time, defaults to last completion',
              default: nil
            },
            override_interval: {
              type: 'integer',
              description: 'Span in seconds, defaults to current time',
              default: nil
            }
          },
          secrets: [
            :sftp_deposit_host, :sftp_deposit_password, :sftp_deposit_user,
            :sftp_ingest_host, :sftp_ingest_password, :sftp_ingest_user
          ],
          workflow_path: '/app/workflows/argo/ingestion/workflow.yaml'
        }
      end
    end

    def validate
    end
  end
end
  
