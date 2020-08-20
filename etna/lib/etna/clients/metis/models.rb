require_relative '../../json_serializable_struct'


module Etna
  module Clients
    class Metis
      class ListFoldersRequest < Struct.new(:project_name, :bucket_name, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({}.update(params))
        end
      end

      class ListFoldersResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def folders
          Folders.new(raw[:folders])
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

      class Folder
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def folder_path
          raw[:folder_path]
        end

        def bucket_name
          raw[:bucket_name]
        end
      end
    end
  end
end