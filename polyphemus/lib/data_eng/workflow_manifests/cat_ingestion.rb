class Polyphemus
    class CatIngestionManifest < Polyphemus::WorkflowManifest
  
      class << self
        def as_json
          {
            name: 'cat',
            schema: {
              regex: 'string',
              root_dir: 'string'
            },
            runtime_params: {
              commit: 'boolean'
            }
            secrets: [ :sftp_host, :sftp_user, :sftp_password, :sftp_port ]
          }
        end
      end
  
      end
  
  
    end
  end
  