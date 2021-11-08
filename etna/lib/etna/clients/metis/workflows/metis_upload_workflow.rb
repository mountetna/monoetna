require 'ostruct'
require 'digest'
require 'fileutils'
require 'tempfile'
require 'securerandom'

$digest_mutx = Mutex.new

module Etna
  module Clients
    class Metis
      class MetisUploadWorkflow < Struct.new(:metis_client, :metis_uid, :project_name, :bucket_name, :max_attempts, keyword_init: true)
        class StreamingUploadError < StandardError
        end


        def initialize(args)
          super({max_attempts: 3, metis_uid: SecureRandom.hex}.update(args))
        end

        def do_upload(source_file_or_upload, dest_path, &block)
          unless source_file_or_upload.is_a?(Upload)
            upload = Upload.new(source_file: source_file_or_upload)
          else
            upload = source_file_or_upload
          end

          dir = ::File.dirname(dest_path)
          metis_client.create_folder(CreateFolderRequest.new(
              project_name: project_name,
              bucket_name: bucket_name,
              folder_path: dir,
          )) unless dir == "."

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
              unless block.nil?
                m = yield [:error, e]
                if m == false
                  raise e
                end
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

        class Upload
          INITIAL_BLOB_SIZE = 2 ** 10
          MAX_BLOB_SIZE = 2 ** 22
          ZERO_HASH = 'd41d8cd98f00b204e9800998ecf8427e'

          attr_accessor :source_file, :next_blob_size, :current_byte_position

          def initialize(source_file: nil, next_blob_size: nil, current_byte_position: nil)
            self.source_file = source_file
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
            bytes = next_blob_bytes
            if bytes.empty?
              return ZERO_HASH
            end

            $digest_mutx.synchronize do
              return Digest::MD5.hexdigest(bytes)
            end
          end

          def next_blob_bytes
            IO.binread(source_file, next_blob_size, current_byte_position)
          end

          def resume_from!(upload_response)
            self.current_byte_position = upload_response.current_byte_position
            self.next_blob_size = upload_response.next_blob_size
          end
        end

        class StreamingIOUpload < Upload
          def initialize(readable_io:, size_hint: 0, **args)
            @readable_io = readable_io
            @size_hint = size_hint
            @read_position = 0
            @last_bytes = ""
            super(**args)
          end

          def file_size
            @size_hint
          end

          def next_blob_bytes
            next_left = current_byte_position
            next_right = current_byte_position + next_blob_size

            if next_right < @read_position
              raise StreamingUploadError.new("Upload needs restart, but source is streaming and ephemeral. #{next_right} #{@read_position} You need to restart the source stream and create a new upload.")
            elsif @read_position < next_left
              # read from the stream and discard until we are positioned for the next read.
              data = @readable_io.read(next_left - @read_position)
              raise StreamingUploadError.new("Unexpected EOF in read stream") if data.nil?

              @read_position += data.bytes.length
            end

            # If we have consumed all requested data, return what we have consumed.
            # If we have requested no data, make sure to provide "" as the result.
            if next_right == @read_position
              return @last_bytes
            end

            if @read_position != next_left
              raise StreamingUploadError.new("Alignment error, source data does not match expected upload resume. #{@read_position} #{next_left} Restart the upload to address.")
            end

            @last_bytes = "".tap do |bytes|
              while @read_position < next_right
                bytes << @readable_io.read(next_right - @read_position).tap do |data|
                  raise StreamingUploadError.new("Unexpected EOF in read stream") if data.nil?
                  @read_position += data.bytes.length
                end
              end
            end
          end
        end
      end
    end
  end
end
