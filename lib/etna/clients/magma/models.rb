require 'ostruct'
require_relative '../../json_serializable_struct'
require_relative '../../multipart_serializable_nested_hash'
require_relative '../../directed_graph'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Magma
      class RetrievalRequest < Struct.new(:model_name, :attribute_names, :record_names, :project_name, :page, :page_size, :order, :filter, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({model_name: 'all', attribute_names: 'all', record_names: []}.update(params))
        end
      end

      class QueryRequest < Struct.new(:query, :project_name, keyword_init: true)
        include JsonSerializableStruct
      end

      class UpdateRequest < Struct.new(:revisions, :project_name, keyword_init: true)
        include JsonSerializableStruct
        include MultipartSerializableNestedHash

        def initialize(**params)
          super({revisions: {}}.update(params))
        end

        def update_revision(model_name, record_name, attrs)
          revision = revisions[model_name] ||= {}
          record = revision[record_name] ||= {}
          record.update(attrs)
        end
      end

      class UpdateModelRequest < Struct.new(:project_name, :actions, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({actions: []}.update(params))
        end

        def add_action(action)
            actions << action
        end
      end

      class AddModelAction < Struct.new(:action_name, :model_name, :parent_model_name, :parent_link_type, :identifier, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({action_name: 'add_model'}.update(args))
        end
      end

      class AddAttributeAction < Struct.new(:action_name, :model_name, :attribute_name, :type, :description, :display_name, :format_hint, :hidden, :index, :link_model_name, :read_only, :restricted, :unique, :validation, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({action_name: 'add_attribute'}.update(args))
        end
      end

      class AddLinkAction < Struct.new(:action_name, :links, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({action_name: 'add_link', links: []}.update(args))
        end
      end

      class AddLinkDefinition < Struct.new(:type, :model_name, :attribute_name, keyword_init: true)
        include JsonSerializableStruct
      end

      class AddProjectAction < Struct.new(:action_name, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({action_name: 'add_project'}.update(args))
        end
      end

      class UpdateAttributeAction < Struct.new(:action_name, :model_name, :attribute_name, :type, :description, :display_name, :format_hint, :hidden, :index, :link_model_name, :read_only, :restricted, :unique, :validation, keyword_init: true)
        include JsonSerializableStruct
        def initialize(**args)
          super({action_name: 'update_attribute'}.update(args))
        end
      end

      class AttributeValidation < Struct.new(:type, :value, :begin, :end, keyword_init: true)
        include JsonSerializableStruct
      end

      class AttributeValidationType < String
        REGEXP = AttributeValidationType.new("Regexp")
        ARRAY = AttributeValidationType.new("Array")
        RANGE = AttributeValidationType.new("Range")
      end

      class RetrievalResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def models
          Models.new(raw['models'])
        end
      end

      class UpdateModelResponse < RetrievalResponse
      end

      class QueryResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def answer
          raw['answer']
        end

        def format
          raw['format']
        end

        def type
          raw['type']
        end
      end

      class UpdateResponse < RetrievalResponse
      end

      class Models
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def model_keys
          raw.keys
        end

        def build_model(model_key)
          Model.new(raw[model_key] ||= {})
        end

        def model(model_key)
          return nil unless raw.include?(model_key)
          Model.new(raw[model_key])
        end

        def to_directed_graph(include_casual_links=false)
          graph = ::DirectedGraph.new

          model_keys.each do |model_name|
            graph.add_connection(model(model_name).template.parent, model_name)

            if include_casual_links
              attributes = model(model_name).template.attributes
              attributes.attribute_keys.each do |attribute_name|
                attribute = attributes.attribute(attribute_name)

                linked_model_name = attribute.link_model_name
                if linked_model_name
                  if attribute.attribute_type == AttributeType::PARENT
                    graph.add_connection(linked_model_name, model_name)
                  elsif attribute.attribute_type == AttributeType::COLLECTION
                    graph.add_connection(model_name, linked_model_name)
                  elsif attribute.attribute_type == AttributeType::CHILD
                    graph.add_connection(model_name, linked_model_name)
                  elsif attribute.attribute_type == AttributeType::LINK
                    graph.add_connection(model_name, linked_model_name)
                  end
                end
              end
            end
          end

          graph
        end
      end

      class Model
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def build_template
          Template.new(raw['template'] ||= {})
        end

        def documents
          Documents.new(raw['documents'])
        end

        def template
          Template.new(raw['template'])
        end

        def count
          raw['count'] || 0
        end
      end

      class Documents
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def document_keys
          raw.keys
        end

        def document(document_key)
          return nil unless raw.include?(document_key)
          raw[document_key]
        end
      end

      class Template
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def name
          raw['name'] || ""
        end

        def name=(val)
          raw['name'] = val.to_s
        end

        def identifier
          raw['identifier'] || ""
        end

        def identifier=(val)
          raw['identifier'] = val.to_s
        end

        def parent
          raw['parent']
        end

        def parent=(val)
          raw['parent'] = val.to_s
        end

        def attributes
          Attributes.new(raw['attributes'])
        end

        def build_attributes
          Attributes.new(raw['attributes'] ||= {})
        end
      end

      class Attributes
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def attribute_keys
          raw.keys
        end

        def attribute(attribute_key)
          return nil unless raw.include?(attribute_key)
          Attribute.new(raw[attribute_key])
        end

        def build_attribute(key)
          Attribute.new(raw[key] ||= {})
        end

        def all
          raw.values.map { |r| Attribute.new(r) }
        end
      end

      class Attribute
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def name
          @raw['name'] || ""
        end

        def name=(val)
          @raw['name'] = val
        end

        def attribute_name
          @raw['attribute_name'] || ""
        end

        def attribute_name=(val)
          @raw['attribute_name'] = val
        end

        def attribute_type
          @raw['attribute_type'] && AttributeType.new(@raw['attribute_type'])
        end

        def attribute_type=(val)
          val = val.to_s if val
          @raw['attribute_type'] = val
        end

        def link_model_name
          @raw['link_model_name']
        end

        def link_model_name=(val)
          @raw['link_model_name'] = val
        end

        def unique
          raw['unique']
        end

        def desc
          raw['desc']
        end

        def desc=(val)
          @raw['desc'] = val
        end

        def display_name
          raw['display_name']
        end

        def display_name=(val)
          raw['display_name'] = val
        end

        def match
          raw['match']
        end

        def restricted
          raw['restricted']
        end

        def format_hint
          raw['format_hint']
        end

        def read_only
          raw['read_only']
        end

        def hidden
          raw['hidden']
        end

        def validation
          raw['validation']
        end

        def options
          raw['options']
        end
      end

      class AttributeType < String
        STRING = AttributeType.new("string")
        DATE_TIME = AttributeType.new("date_time")
        BOOLEAN = AttributeType.new("boolean")
        CHILD = AttributeType.new("child")
        COLLECTION = AttributeType.new("collection")
        FILE = AttributeType.new("file")
        FLOAT = AttributeType.new("float")
        IDENTIFIER = AttributeType.new("identifier")
        IMAGE = AttributeType.new("image")
        INTEGER = AttributeType.new("integer")
        LINK = AttributeType.new("link")
        MATCH = AttributeType.new("match")
        MATRIX = AttributeType.new("matrix")
        PARENT = AttributeType.new("parent")
        TABLE = AttributeType.new("table")
      end

      class ParentLinkType < String
        CHILD = ParentLinkType.new("child")
        COLLECTION = ParentLinkType.new("collection")
        TABLE = ParentLinkType.new("table")
      end
    end
  end
end