require 'struct'

module Etna
  module Clients
    module Magma
      class RetrievalRequest < Struct.new(:model_name, :attribute_names, :record_names, :project_name, keyword_init: true)
      end

      class RetrievalResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = {'models': {}}.update(raw)
        end

        def models
          Models.new(raw['models'])
        end
      end

      class Models
        attr_reader :raw

        def initialize(raw = {})
          @raw = {}.update(raw)
        end

        def model_keys
          raw.keys
        end

        def model(model_key)
          Model.new(raw[model_key])
        end
      end

      class Model
        attr_reader :raw

        def initialize(raw = {})
          @raw = {'documents': {}}.update(raw)
        end

        def documents
          Documents.new(raw['documents'])
        end
      end

      class Documents
        attr_reader :raw

        def initialize(raw = {})
          @raw = {}.update(raw)
        end

        def document_keys
          raw.keys
        end

        def document(document_key)
          raw[document_key]
        end
      end

      class Template
        attr_reader :raw

        def initialize(raw = {})
          @raw = {'attributes': {}}.update(raw)
        end

        def name
          raw['name'] || ""
        end

        def identifier
          raw['identifier'] || ""
        end

        def parent
          raw['parent']
        end

        def attributes
          Attributes.new(raw['attributes'])
        end
      end

      class Attributes
        attr_reader :raw

        def initialize(raw = {})
          @raw = {}.update(raw)
        end

        def attribute_keys
          raw.keys
        end

        def attribute(attribute_key)
          Attribute.new(raw[attribute_key])
        end
      end

      class Attribute
        attr_reader :raw

        def initialize(raw = {})
          @raw = {}.update(raw)
        end

        def name
          @raw['name'] || ""
        end

        def attribute_name
          @raw['attribute_name'] || ""
        end

        def type
          @raw['type'] && Type.new(@raw['type'])
        end
      end

      class Type < String
        STRING = Type.new("string")
        DATE_TIME = Type.new("date_time")
        BOOLEAN = Type.new("boolean")
        CHILD = Type.new("child")
        COLLECTION = Type.new("collection")
        FILE = Type.new("file")
        FLOAT = Type.new("float")
        IDENTIFIER = Type.new("identifier")
        IMAGE = Type.new("image")
        INTEGER = Type.new("integer")
        LINK = Type.new("link")
        MATCH = Type.new("match")
        MATRIX = Type.new("matrix")
        PARENT = Type.new("parent")
        TABLE = Type.new("table")
      end
    end
  end
end