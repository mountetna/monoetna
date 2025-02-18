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
              root_dir: { type: 'string' },
              file_regex: { type: 'string' },
              sftp_root_dir: { type: 'string' },
              path_to_write_files: { type: 'string' },
              bucket_name: { type: 'string' },
              metis_root_path: { type: 'string' },
              c4_root_path: { type: 'string' }
            }
          },
          runtime_params: {
            commit: 'boolean',
            initial_start_scan_time: 'integer', # unix timestamp
          },
          secrets: [:sftp_host, :sftp_user, :sftp_password, :sftp_port],
          workflow_path: '/app/workflows/argo/ingestion/workflow.yaml'
        }
      end
    end

    def validate
    end
  end
end
  
