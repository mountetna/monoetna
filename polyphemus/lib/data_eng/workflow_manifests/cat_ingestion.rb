class Polyphemus
    class CatIngestionManifest < Polyphemus::WorkflowManifest
  
      class << self
        def as_json
          {
            name: 'cat',
              schema: {
                regex: 'string',
                root_dir: 'string',
                file_regex: 'string',
                sftp_root_dir: 'string',
                interval: 'integer', # in seconds 
                initial_start_scan_time: 'integer', # unix timestamp
                path_to_write_files: 'string',
                bucket_name: 'string',
                metis_root_path: 'string',
                c4_root_path: 'string'
                },
            runtime_params: {
              commit: 'boolean'
            },
            secrets: [:sftp_host, :sftp_user, :sftp_password, :sftp_port]
          }
        end
      end
  
    end
  end
end
  