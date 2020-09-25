require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'
require 'securerandom'

module Etna
  module Clients
    class Metis
      class MetisUploadWorkflow < Struct.new(:metis_client, :metis_uid, :project_name, :bucket_name, :max_attempts, keyword_init: true)

        def initialize(args)
          super({max_attempts: 3, metis_uid: SecureRandom.hex}.update(args))
        end

        def do_upload(source_file, dest_path, &block)
          upload = Upload.new(source_file: source_file)

          dir = ::File.dirname(dest_path)
          metis_client.create_folder(CreateFolderRequest.new(
              project_name: project_name,
              bucket_name: bucket_name,
              folder_path: dir,
          ))

          authorize_response = metis_client.authorize_upload(AuthorizeUploadRequest.new(
              project_name: project_name,
              bucket_name: bucket_name,
              file_path: dest_path,
          ))

          upload_parts(upload, authorize_response.upload_path, &block)
        end

        private

        def upload_parts(upload, upload_path, attempt_number = 1, reset = false, &block)
          if attempt_number > max_attempts
            raise "Upload failed after (#{attempt_number}) attempts."
          end

          start = Time.now
          unsent_zero_byte_file = upload.current_byte_position == 0

          upload.resume_from!(metis_client.upload_start(UploadStartRequest.new(
              upload_path: upload_path,
              file_size: upload.file_size,
              next_blob_size: upload.next_blob_size,
              next_blob_hash: upload.next_blob_hash,
              metis_uid: metis_uid,
              reset: reset,
          )))

          until upload.complete? && !unsent_zero_byte_file
            begin
              blob_bytes = upload.next_blob_bytes
              byte_position = upload.current_byte_position

              upload.advance_position!
              metis_client.upload_blob(UploadBlobRequest.new(
                  upload_path: upload_path,
                  next_blob_size: upload.next_blob_size,
                  next_blob_hash: upload.next_blob_hash,
                  blob_data: StringIO.new(blob_bytes),
                  metis_uid: metis_uid,
                  current_byte_position: byte_position,
              ))

              unsent_zero_byte_file = false
            rescue Etna::Error => e
              m = yield [:error, e] unless block.nil?
              if m == false
                raise e
              end

              if e.status == 422
                return upload_parts(upload, upload_path, attempt_number + 1, true, &block)
              elsif e.status >= 500
                return upload_parts(upload, upload_path, attempt_number + 1, &block)
              end

              raise e
            end

            yield [
                :progress,
                upload.file_size == 0 ? 1.0 : upload.current_byte_position.to_f / upload.file_size,
                (upload.current_byte_position / (Time.now - start).to_f).round(2),
            ] unless block.nil?
          end
        end

        class Upload < Struct.new(:source_file, :next_blob_size, :current_byte_position, keyword_init: true)
          INITIAL_BLOB_SIZE = 2 ** 10
          MAX_BLOB_SIZE = 2 ** 22
          ZERO_HASH = 'd41d8cd98f00b204e9800998ecf8427e'

          def initialize(**args)
            super
            self.next_blob_size = [file_size, INITIAL_BLOB_SIZE].min
            self.current_byte_position = 0
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
            Digest::MD5.hexdigest(next_blob_bytes)
          end

          def next_blob_bytes
            IO.binread(source_file, next_blob_size, current_byte_position)
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
