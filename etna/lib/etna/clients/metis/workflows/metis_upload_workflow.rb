require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'

module Etna
  module Clients
    class Metis
      class MetisUploadWorkflow < Struct.new(:metis_client, :project_name, :bucket_name, :max_attempts, keyword_init: true)

        def initialize(args)
          super({max_attempts: 3}.update(args))
        end

        def do_upload(source_file, dest_project, dest_bucket, dest_path, &block)
          upload = Upload.new(source_file: source_file)

          authorize_response = metis_client.authorize_upload(AuthorizeUploadRequest.new(
              project_name: dest_project,
              bucket_name: dest_bucket,
              file_path: dest_path,
          ))

          upload.resume_from!(metis_client.upload_start(UploadStartRequest.new(
              file_size: upload.file_size,
              next_blob_size: upload.next_blob_size,
              next_blob_hash: upload.next_blob_hash,
              upload_path: authorize_response.upload_path,
          )))

          max_attempts.times do
            if upload.complete?
              metis_client.upload_blob(UploadBlobRequest.new(
                  next_blob_size: 0,
                  next_blob_hash: ZERO_HASH,
                  upload_path: authorize_response.upload_path,
                  blob_data: '',
                  current_byte_position: upload.current_byte_position
              ))

              yield [1, 0]
              return
            end
          end
        end

        class Upload < Struct.new(:source_file, :next_blob_size, :current_byte_position, keyword_init: true)
          INITIAL_BLOB_SIZE = 2 ** 10
          MAX_BLOB_SIZE = 2 ** 22
          ZERO_HASH = 'd41d8cd98f00b204e9800998ecf8427e'

          def initialize(file_size:, **args)
            super({next_blob_size: [file_size, INITIAL_BLOB_SIZE].min, current_byte_position: 0}.update(args))
          end

          def file_size
            ::File.size(source_file)
          end

          def advance_position!
            self.current_byte_position = self.current_byte_position + self.next_blob_size
            self.next_blob_size = [
                                  MAX_BLOB_SIZE,
                                  # in fact we should stop when we hit the end of the file
                                  file_size - current_byte_position
                              ].min
          end

          def complete?
            current_byte_position >= file_size
          end

          def next_blob_hash
            Digest::MD5.hexdigest(current_bytes)
          end

          def current_bytes
            IO.binread(file_path, next_blob_size, current_byte_position)
          end

          def resume_from!(upload_response)
            self.current_byte_position = upload_response.current_byte_position
            self.next_blob_size = upload_response.next_blob_size
          end
        end
      end
    end
  end
end
