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

      class FolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, keyword_init: true)
      end

      class ListFolderRequest < FolderRequest
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

      class ListFolderByIdRequest < Struct.new(:project_name, :bucket_name, :folder_id, keyword_init: true)
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

      class TouchFolderRequest < Struct.new(:project_name, :bucket_name, :folder_path, keyword_init: true)
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

      class TouchFileRequest < Struct.new(:project_name, :bucket_name, :file_path, keyword_init: true)
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

      class CreateFolderRequest < FolderRequest
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

      class DeleteFolderRequest < FolderRequest
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

      class TailBucketRequest < Struct.new(:project_name, :bucket_name, :type, :folder_id, :batch_start, :batch_end, keyword_init: true)
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

      class TailNode < Etna::Clients::Response
        def initialize(tail, raw = {})
          super(raw)
          @tail = tail
        end

        def parent?
          @raw[:type] == 'parent'
        end

        def folder?
          @raw[:type] == 'folder'
        end

        def file?
          @raw[:type] == 'file'
        end

        params :id, :node_name, :parent_id, :file_hash, :updated_at

        def as_file
          {
            file_name: node_name,
            file_hash: file_hash,
            updated_at: updated_at,
            file_path: file_path,
            folder_id: parent_id,
            project_name: @tail.project_name,
            bucket_name: @tail.bucket_name,
          }
        end

        def parent_path
          return nil if !parent_id

          if !@tail.paths[parent_id]
            parent = @tail.parents[parent_id]
            @tail.paths[parent_id] = parent.file_path
          end

          return @tail.paths[parent_id]
        end

        def file_path
          [ parent_path, node_name ].compact.join('/')
        end
      end

      class TailResponse < Etna::Clients::Response
        attr_reader :project_name, :bucket_name, :paths, :parents


        def initialize(project_name, bucket_name, raw = {})
          super(raw)

          @paths = {}

          @parents = {}

          @project_name = project_name
          @bucket_name = bucket_name

          @nodes = []
          raw.each do |raw_node|
            node = TailNode.new(self, raw_node)
            if node.parent?
              @parents[node.id] = node
            else
              @nodes.push(node)
            end
          end
        end

        def folders
          @nodes.filter_map do |node|
            Folder.new(node.as_folder) if node.folder?
          end
        end

        def files
          @nodes.filter_map do |node|
            File.new(node.as_file) if node.file?
          end
        end
      end

      class FoldersResponse < Etna::Clients::Response
        def folders
          Folders.new(raw[:folders])
        end
      end

      class FilesResponse < Etna::Clients::Response
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

      class Files < Etna::Clients::Response
        def all
          raw.map { |file| File.new(file) }
        end
      end

      class Folders < Etna::Clients::Response
        def all
          raw.map { |folder| Folder.new(folder) }
        end
      end

      class File < Etna::Clients::Response
        def with_containing_folder(folder)
          folder_path = folder.is_a?(Folder) ? folder.folder_path : folder
          File.new({}.update(self.raw).update({
            file_path: ::File.join(folder_path, self.file_name)
          }))
        end

        params :file_path, :project_name, :bucket_name, :file_name, :size, :file_hash, :folder_id

        def download_path
          raw[:download_url].nil? ?
              "/#{project_name}/download/#{bucket_name}/#{file_path}" :
              raw[:download_url].sub(%r!^https://[^/]*?/!, '/')
        end

        def download_url
          raw[:download_url] || ''
        end

        def updated_at
          time = raw[:updated_at]
          time.nil? ? nil : Time.parse(time)
        end

        def as_magma_file_attribute
          {
            path: as_metis_url,
            original_filename: file_name || ::File.basename(file_path)
          }
        end

        def as_metis_url
          "metis://#{project_name}/#{bucket_name}/#{file_path}"
        end
      end

      class Folder < Etna::Clients::Response
        params :folder_path, :folder_name, :bucket_name, :project_name, :id

        def updated_at
          time = raw[:updated_at]
          time.nil? ? nil : Time.parse(time)
        end
      end

      class AuthorizeUploadRequest < Struct.new(:project_name, :bucket_name, :file_path, keyword_init: true)
        include JsonSerializableStruct
      end

      class AuthorizeDownloadRequest < Struct.new(:project_name, :bucket_name, :file_path, keyword_init: true)
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

      class DownloadResponse < Etna::Clients::Response
        def download_url
          raw['download_url']
        end
      end

      class UploadResponse < Etna::Clients::Response
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

      class GetFileCountByProjectRequest < Struct.new(:project_names, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end
      end

      class GetByteCountByProjectRequest < Struct.new(:project_names, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end
      end
    end
  end
end
