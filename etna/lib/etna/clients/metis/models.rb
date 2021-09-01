require 'ostruct'
require_relative '../../json_serializable_struct'
require_relative '../base_client'

module Etna
  module Clients
    class Metis < Etna::Clients::BaseClient
      class ListFoldersRequest < Struct.new(:project_name, :bucket_name, :offset, :limit, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class RenameFolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, :new_bucket_name, :new_folder_path, :create_parent, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({create_parent: false}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class RenameFileRequest < Struct.new(:project_name, :bucket_name, :file_path, :new_bucket_name, :new_file_path, :create_parent, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({create_parent: false}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class ListFolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class CreateFolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class DeleteFolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class DeleteFileRequest < Struct.new(:project_name, :bucket_name, :file_path, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end

        def to_h
          # The :project_name comes in from Polyphemus as a symbol value,
          #   we need to make sure it's a string because it's going
          #   in the URL.
          super().compact.transform_values(&:to_s)
        end
      end

      class FindRequest < Struct.new(:project_name, :bucket_name, :limit, :offset, :params, :hide_paths, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({params: [], hide_paths: false}.update(args))
        end

        def add_param(param)
          params << param
        end

        def to_h
          # The nested :params values don't get converted correctly with transform_values, so it's
          #   easier to do from a JSON string
          JSON.parse(to_json, :symbolize_names => true)
        end

        def clone
          FindRequest.new(
            project_name: self.project_name,
            bucket_name: self.bucket_name,
            limit: self.limit,
            offset: self.offset,
            params: self.params.dup,
            hide_paths: self.hide_paths
          )
        end
      end

      class FindParam < Struct.new(:attribute, :predicate, :value, :type, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({}.update(args))
        end
      end

      class CopyFilesRequest < Struct.new(:project_name, :revisions, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({revisions: []}.update(args))
        end

        def add_revision(revision)
          revisions << revision
        end

        def to_h
          # The nested :revisions values don't get converted correctly with transform_values, so it's
          #   easier to do from a JSON string
          JSON.parse(to_json, :symbolize_names => true)
        end
      end

      class CopyRevision < Struct.new(:source, :dest, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({}.update(args))
        end
      end

      class FoldersResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def folders
          Folders.new(raw[:folders])
        end
      end

      class FilesResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def files
          Files.new(raw[:files])
        end
      end

      class FoldersAndFilesResponse < FoldersResponse
        def files
          Files.new(raw[:files] || [])
        end

        def folders
          Folders.new(raw[:folders] || [])
        end
      end

      class Files
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def all
          raw.map { |file| File.new(file) }
        end
      end

      class Folders
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def all
          raw.map { |folder| Folder.new(folder) }
        end
      end

      class File
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def file_path
          raw[:file_path]
        end

        def project_name
          raw[:project_name]
        end

        def bucket_name
          raw[:bucket_name]
        end

        def download_path
          raw[:download_url].nil? ?
              "/#{project_name}/download/#{bucket_name}/#{file_path}" :
              raw[:download_url].sub(%r!^https://[^/]*?/!, '/')
        end

        def download_url
          raw[:download_url] || ''
        end

        def file_name
          raw[:file_name]
        end

        def updated_at
          time = raw[:updated_at]
          time.nil? ? nil : Time.parse(time)
        end

        def size
          raw[:size]
        end

        def file_hash
          raw[:file_hash]
        end

        def folder_id
          raw[:folder_id]
        end
      end

      class Folder
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def folder_path
          raw[:folder_path]
        end

        def folder_name
          raw[:folder_name]
        end

        def bucket_name
          raw[:bucket_name]
        end

        def updated_at
          time = raw[:updated_at]
          time.nil? ? nil : Time.parse(time)
        end

        def id
          raw[:id]
        end
      end

      class AuthorizeUploadRequest < Struct.new(:project_name, :bucket_name, :file_path, keyword_init: true)
        include JsonSerializableStruct
      end

      class UploadStartRequest < Struct.new(:file_size, :action, :metis_uid, :next_blob_size, :upload_path, :next_blob_hash, :reset, keyword_init: true)
        include JsonSerializableStruct

        def initialize(args)
          super({ action: UploadAction::START }.update(args))
        end
      end

      class UploadBlobRequest < Struct.new(:file_size, :action, :metis_uid, :blob_data, :upload_path, :next_blob_size, :next_blob_hash, :current_byte_position, keyword_init: true)
        include MultipartSerializableNestedHash

        def initialize(args)
          super({ action: UploadAction::BLOB }.update(args))
        end

        def encode_multipart_content(base_key = '')
          self.class.encode_multipart_content(to_h, base_key)
        end
      end

      class UploadResponse
        attr_reader :raw
        def initialize(raw = {})
          @raw = raw
        end

        def current_byte_position
          raw['current_byte_position'].to_i
        end

        def url
          raw['url'] || ''
        end

        def next_blob_size
          raw['next_blob_size'].to_i
        end

        def upload_path
          url.sub(%r!^https://[^/]*?/!, '/')
        end
      end

      class UploadAction < String
        START = UploadAction.new("start")
        BLOB = UploadAction.new("blob")
      end
    end
  end
end