require 'ostruct'
require_relative '../../json_serializable_struct'
require_relative '../../multipart_serializable_nested_hash'
require_relative '../../directed_graph'
require_relative '../enum'
require_relative '../base_client'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Magma < Etna::Clients::BaseClient
      class RetrievalRequest < Struct.new(:model_name, :attribute_names, :record_names, :project_name, :page, :page_size, :order, :filter, :hide_templates, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({model_name: 'all', attribute_names: 'all', record_names: []}.update(params))
        end
      end

      class QueryRequest < Struct.new(:query, :project_name, :order, :page, :page_size, keyword_init: true)
        include JsonSerializableStruct
      end

      class UpdateRequest < Struct.new(:revisions, :project_name, :dry_run, keyword_init: true)
        include JsonSerializableStruct
        include MultipartSerializableNestedHash

        def initialize(**params)
          super({revisions: {}, dry_run: false}.update(params))
        end

        def update_revision(model_name, record_name, attrs)
          revision = revisions[model_name] ||= {}
          record = revision[record_name] ||= {}
          record.update(attrs)
        end

        def append_table(parent_model_name, parent_record_name, model_name, attrs, attribute_name = model_name)
          parent_revision = update_revision(parent_model_name, parent_record_name, {})
          table = parent_revision[attribute_name] ||= []
          id = "::#{model_name}#{(revisions[model_name] || {}).length + 1}"
          table << id
          update_revision(model_name, id, attrs)
          id
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

      class AddModelAction < Struct.new(:action_name, :model_name, :parent_model_name, :parent_link_type, :identifier, :date_shift_root, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'add_model', date_shift_root: false}.update(args))
        end
      end

      class SetDateShiftRootAction < Struct.new(:action_name, :model_name, :date_shift_root, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'set_date_shift_root'}.update(args))
        end
      end

      class AddAttributeAction < Struct.new(:action_name, :model_name, :attribute_name, :type, :description, :display_name, :format_hint, :hidden, :index, :link_model_name, :read_only, :attribute_group, :restricted, :unique, :validation, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'add_attribute'}.update(args))
        end

        def attribute_type=(val)
          self.type = val
        end

        def attribute_type
          self.type
        end

        def desc=(val)
          self.description = val
        end

        def desc
          self.description
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

      class AddProjectAction < Struct.new(:action_name, :no_metis_bucket, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'add_project'}.update(args))
        end
      end

      class UpdateAttributeAction < Struct.new(:action_name, :model_name, :attribute_name, :description, :display_name, :format_hint, :hidden, :index, :link_model_name, :read_only, :attribute_group, :restricted, :unique, :validation, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'update_attribute'}.update(args))
        end

        def as_json
          super(keep_nils: true)
        end
      end

      class RenameAttributeAction < Struct.new(:action_name, :model_name, :attribute_name, :new_attribute_name, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**args)
          super({action_name: 'rename_attribute'}.update(args))
        end
      end

      class AttributeValidation < Struct.new(:type, :value, :begin, :end, keyword_init: true)
        include JsonSerializableStruct
      end

      class AttributeValidationType < String
        include Enum
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

      class Project
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def models
          Models.new(raw['models'])
        end
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

        def all
          raw.values.map { |r| Model.new(r) }
        end

        def +(other)
          raw_update = {}
          raw_update[other.name] = other.raw
          Models.new({}.update(raw).update(raw_update))
        end

        # Can look up reciprocal links by many means.  At minimum, a source model or model name must be provided,
        # and either a link_attribute_name, attribute, or link_model.
        def find_reciprocal(
            model_name: nil,
            model: self.model(model_name),
            link_attribute_name: nil,
            attribute: model&.template&.attributes&.attribute(link_attribute_name),
            link_model: self.model(attribute&.link_model_name)
        )
          return nil if model.nil? || model.name.nil?
          link_model&.template&.attributes&.all&.find { |a| a.link_model_name == model.name }
        end

        def to_directed_graph(include_casual_links = false)
          graph = ::DirectedGraph.new

          model_keys.sort.each do |model_name|
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

        def name
          @raw.dig('template', 'name')
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

        def +(other)
          Documents.new({}.update(raw).update(other.raw))
        end

        def document(document_key)
          if document_key.is_a?(String)
            raw[document_key]
          else
            raw[document_key&.to_s]
          end
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

        def version
          raw['version'] || 0
        end

        def version=(val)
          raw['version'] = val
        end

        def identifier=(val)
          raw['identifier'] = val
        end

        def parent
          raw['parent']
        end

        def parent=(val)
          raw['parent'] = val
        end

        def attributes
          Attributes.new(raw['attributes'] ||= {})
        end

        def build_attributes
          Attributes.new(raw['attributes'] ||= {})
        end

        def dictionary
          Dictionary.new(raw['dictionary'] ||= {})
        end

        def build_dictionary
          Dictionary.new(raw['dictionary'] ||= {})
        end

        def all_linked_model_names
          models = [ self.parent, ] + build_attributes.all.map { |v| v.link_model_name }
          models.select { |m| !m.nil? }.uniq
        end
      end

      class Dictionary
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def dictionary_keys
          raw.keys
        end

        def dictionary_model
          raw['dictionary_model']
        end

        def dictionary_model=(val)
          @raw['dictionary_model'] = val
        end

        def model_name
          raw['model_name']
        end

        def model_name=(val)
          @raw['model_name'] = val
        end

        def attributes
          raw['attributes']
        end

        def attributes=(val)
          @raw['attributes'] = val
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

        def is_edited?(other)
          # Don't just override == in case need to do a full comparison.
          editable_attribute_names = Attribute::EDITABLE_ATTRIBUTE_ATTRIBUTES.map(&:to_s)
          
          self_editable = raw.slice(*editable_attribute_names)
          other_editable = other.raw.slice(*editable_attribute_names)

          self_editable != other_editable
        end

        # Sets certain attribute fields which are implicit, even when not set, to match server behavior.
        def set_field_defaults!
          @raw.replace({
              'hidden' => false,
              'read_only' => false,
              'restricted' => false,
              'validation' => nil,
          }.merge(@raw))
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

        def is_project_name_reference?(model_name)
          return true if model_name == 'project' && attribute_type == AttributeType::IDENTIFIER
          return true if link_model_name == 'project'
          false
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

        def unique=(val)
          raw['unique'] = val
        end

        def description
          raw['description']
        end

        def description=(val)
          @raw['description'] = val
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

        def restricted=(val)
          raw['restricted'] = val
        end

        def format_hint
          raw['format_hint']
        end

        def format_hint=(val)
          raw['format_hint'] = val
        end

        def read_only
          raw['read_only']
        end

        def read_only=(val)
          raw['read_only'] = val
        end

        def attribute_group
          raw['attribute_group']
        end

        def attribute_group=(val)
          raw['attribute_group'] = val
        end

        def hidden
          raw['hidden']
        end

        def hidden=(val)
          raw['hidden'] = val
        end

        def validation
          raw['validation']
        end

        def validation=(val)
          raw['validation'] = val
        end

        def options
          raw['options']
        end

        COPYABLE_ATTRIBUTE_ATTRIBUTES = [
            :attribute_name, :attribute_type, :display_name, :format_hint,
            :hidden, :link_model_name, :read_only, :attribute_group, :unique, :validation,
            :restricted, :description
        ]

        EDITABLE_ATTRIBUTE_ATTRIBUTES = UpdateAttributeAction.members & COPYABLE_ATTRIBUTE_ATTRIBUTES

        def self.copy(source, dest, attributes: COPYABLE_ATTRIBUTE_ATTRIBUTES)
          attributes.each do |attr_name|
            next unless dest.respond_to?(:"#{attr_name}=")
            source_v = source.send(attr_name)
            dest.send(:"#{attr_name}=", source_v)
          end
        end
      end

      class AttributeType < String
        include Enum
        STRING = AttributeType.new("string")
        DATE_TIME = AttributeType.new("date_time")
        BOOLEAN = AttributeType.new("boolean")
        CHILD = AttributeType.new("child")
        COLLECTION = AttributeType.new("collection")
        FILE = AttributeType.new("file")
        FILE_COLLECTION = AttributeType.new("file_collection")
        FLOAT = AttributeType.new("float")
        IDENTIFIER = AttributeType.new("identifier")
        IMAGE = AttributeType.new("image")
        INTEGER = AttributeType.new("integer")
        LINK = AttributeType.new("link")
        MATCH = AttributeType.new("match")
        MATRIX = AttributeType.new("matrix")
        PARENT = AttributeType.new("parent")
        TABLE = AttributeType.new("table")
        SHIFTED_DATE_TIME = AttributeType.new("shifted_date_time")
      end

      class ParentLinkType < String
        include Enum
        CHILD = ParentLinkType.new("child")
        COLLECTION = ParentLinkType.new("collection")
        TABLE = ParentLinkType.new("table")
      end
    end
  end
end
