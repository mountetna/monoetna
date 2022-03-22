require 'yaml'

module Etna
  class Cwl
    FIELD_LOADERS = {}

    def initialize(attributes)
      @attributes = attributes
    end

    def self.loader
      Etna::Cwl::RecordLoader.new(self)
    end

    def self.as_json(obj)
      if obj.is_a?(Cwl)
        as_json(obj.instance_variable_get(:@attributes))
      elsif obj.is_a?(Hash)
        {}.tap do |result|
          obj.each do |k, v|
            result[k] = as_json(v)
          end
        end
      elsif obj.is_a?(Array)
        obj.map { |v| as_json(v) }
      else
        obj
      end
    end

    def as_json
      self.class.as_json(@attributes)
    end

    class Loader
      def load(val)
        raise "Unimplemented"
      end

      def optional
        OptionalLoader.new(self)
      end

      def map(&block)
        FunctorMapLoader.new(self, &block)
      end

      def as_mapped_array(id_key = nil, value_key = nil)
        MapLoader.new(self.as_array, id_key, value_key)
      end

      def or(*alternatives)
        UnionLoader.new(self, *alternatives)
      end

      def as_array
        ArrayLoader.new(self)
      end
    end

    class AnyLoader < Loader
      def load(val)
        val
      end

      ANY = AnyLoader.new
    end

    class PrimitiveLoader < Loader
      def initialize(name, type)
        @name = name
        @type = type
      end

      def load(val)
        unless val.is_a?(@type)
          raise "Unexpected val #{val.inspect} for #{@name} type"
        end

        val
      end

      def name
        @name
      end

      def self.find_primitive_type_loader(type_name)
        constants.each do |c|
          c = const_get(c)
          if c.is_a?(Loader)
            return c if c.name == type_name
          end
        end
      end

      STRING = PrimitiveLoader.new('string', String)
      INT = PrimitiveLoader.new('int', Integer)
      LONG = PrimitiveLoader.new('long', Integer)
      FLOAT = PrimitiveLoader.new('float', Float)
      DOUBLE = PrimitiveLoader.new('double', Float)
      NULL = PrimitiveLoader.new('null', NilClass)

      class BooleanLoader < Loader
        def name
          'boolean'
        end

        def load(val)
          raise "Invalid value #{val.inspect} for boolean" unless val.instance_of?(TrueClass) || val.instance_of?(FalseClass)
          val
        end
      end

      BOOLEAN = BooleanLoader.new
    end

    class SourceLoader < Loader
      # Resolves a string of the forms "a-primary-identifier" or "step_name/output" into
      #   [:primary_inputs, "a-primary-identifier"] or
      #   ["step_name", "output"] respectively
      def load(val)
        parts = []

        if val.is_a?(Symbol)
          val = val.to_s
        end

        if val.is_a?(Array)
          parts = PrimitiveLoader::STRING.as_array.load(val)
        elsif val.is_a?(String)
          parts = val.split('/', max = 2)
        end

        if parts.length == 1
          return [:primary_inputs, parts[0]]
        elsif parts.length == 2
          return parts
        end

        raise "Unexpected value for source #{val.inspect}"
      end
    end

    class StrictMapLoader < Loader
      def initialize(items, keys)
        @items = items
        @keys = keys
      end

      def load(val)
        if val.is_a?(Hash)
          val.map do |k, v|
            [@keys.load(k), @items.load(v)]
          end.to_h
        else
          raise "Unexpected val #{val.inspect} for hash"
        end
      end
    end

    class MapLoader < Loader
      def initialize(items, idKey = nil, valueKey = nil)
        @items = items
        @idKey = idKey
        @valueKey = valueKey
      end

      def load(val)
        if val.is_a?(Hash)
          val = [].tap do |result|
            errors = {}
            val.keys.sort.each do |k|
              begin
                v = val[k]
                if v.is_a?(Hash)
                  v[@idKey] = k
                else
                  v = {@idKey => k, @valueKey => v}
                end

                result << v
              rescue => e
                errors[k] = e.to_s
              end
            end

            unless errors.empty?
              raise errors.map { |k, v| "#{k}: #{v}" }.join("\n")
            end
          end
        end

        @items.load(val)
      end
    end

    class ArrayLoader < Loader
      def initialize(items)
        @items = items
      end

      def load(val)
        unless val.is_a?(Array)
          raise "Unexpected val #{val.inspect} for array"
        end

        [].tap do |result|
          errors = []
          val.each do |item|
            begin
              loaded = Cwl.load_item(item, UnionLoader.new(self, @items))
              if loaded.is_a?(Array)
                result.push(*loaded)
              else
                result << loaded
              end
            rescue => e
              errors << e.to_s
            end
          end

          unless errors.empty?
            raise errors.join("\n")
          end
        end
      end
    end

    class EnumLoader < Loader
      def initialize(*options)
        @options = options
      end

      def load(val)
        if @options.include?(val)
          return val
        end

        raise "Value #{val.inspect} does not belong to one of (#{@options.join(', ')})"
      end

      PRIMITIVE_TYPE = EnumLoader.new("null", "boolean", "int", "long", "float", "double", "string")
      NOMINAL_TYPE = EnumLoader.new("File")
    end

    class OptionalLoader < Loader
      def initialize(inner_loader)
        @inner_loader = inner_loader
      end

      def load(val)
        if val.nil?
          return nil
        end

        @inner_loader.load(val)
      end
    end

    class RecordLoader < Loader
      def initialize(klass, field_loaders = nil)
        @klass = klass
        @field_loaders = field_loaders
      end

      def field_loaders
        @field_loaders || @klass::FIELD_LOADERS
      end

      def load(val)
        unless val.is_a?(Hash)
          raise "Unexpected value #{val.inspect} for type #{@klass.name}"
        end

        errors = {}
        @klass.new({}.tap do |result|
          field_loaders.each do |field_sym, loader|
            field_str = field_sym.to_s
            begin
              result[field_str] = loader.load(val[field_str])
            rescue => e
              errors[field_str] = e.to_s
            end
          end

          unless errors.empty?
            raise errors.map { |k, e| "#{k}: #{e}" }.join(',')
          end
        end)
      end
    end

    class NeverLoader < Loader
      def load(val)
        raise "This feature is not supported"
      end

      UNSUPPORTED = NeverLoader.new.optional
    end

    class UnionLoader < Loader
      def initialize(*alternatives)
        @alternatives = alternatives
      end

      def load(val)
        errors = []
        @alternatives.each do |loader|
          begin
            return loader.load(val)
          rescue => e
            errors << e.to_s
          end
        end

        raise errors.join(", ")
      end
    end

    class ArrayType < Cwl
      class InnerLoader < Loader
        def load(val)
          RecordLoader.new(ArrayType, {
              type: EnumLoader.new("array"),
              items: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
          }).load(val)
        end
      end

      def type_loader
        loader = RecordType::Field.type_loader(@attributes['items'])
        return nil if loader.nil?
        @type_loader ||= ArrayLoader.new(loader)
      end
    end

    class EnumType < Cwl
      class InnerLoader < Loader
        def load(val)
          RecordLoader.new(EnumType, {
              type: EnumLoader.new("enum"),
              symbols: PrimitiveLoader::STRING.as_array,
          }).load(val)
        end
      end

      def type_loader
        @type_loader ||= EnumLoader.new(*@attributes['symbols'])
      end
    end

    class RecordType < Cwl
      class RecordTypeLoader
        def load(val)
          RecordLoader.new(RecordType, {
              type: EnumLoader.new("record"),
              fields: Field::FieldLoader.new.as_mapped_array('name', 'type')
          }).load(val)
        end
      end

      class Record
        def self.new(h)
          h
        end
      end

      def type_loader
        @type_loader ||= begin
          record_class = Class.new(Record)
          RecordLoader.new(record_class, @attributes['fields'].map do |field|
            loader = Field.type_loader(field.type)
            return nil if loader.nil?
            [field.name.to_sym, loader]
          end.to_h)
        end
      end

      class Field < Cwl
        class FieldLoader < Loader
          def load(val)
            RecordLoader.new(Field, {
                name: PrimitiveLoader::STRING,
                type: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
                doc: PrimitiveLoader::STRING.optional,
            }).load(val)
          end
        end

        def name
          @attributes['name']
        end

        def type
          @attributes['type']
        end

        def type_loader
          self.class.type_loader(self.type)
        end

        def self.type_loader(type)
          case type
          when Array
            type_loaders = type.map { |t| Field.type_loader(t) }
            return nil if type_loaders.any?(&:nil?)
            UnionLoader.new(*type_loaders)
          when EnumType
            type.type_loader
          when RecordType
            type.type_loader
          when ArrayType
            type.type_loader
          when String
            PrimitiveLoader.find_primitive_type_loader(type)
          else
            raise "Could not determine loader for type #{type.inspect}"
          end
        end
      end
    end

    class FunctorMapLoader < Loader
      def initialize(inner, &block)
        @block = block
        @inner = inner
      end

      def load(val)
        @block.call(@inner.load(val))
      end
    end

    # Prepares a unique set of structured nominal types for an inner
    # loading of types
    class TypedDSLLoader < Loader
      def initialize(inner)
        @inner = inner
      end

      REGEX = /^([^\[?]+)(\[\])?(\?)?$/

      def resolve(val)
        m = REGEX.match(val)

        unless m.nil?
          type = m[1]
          unless m[2].nil?
            type = {'type' => 'array', 'items' => type}
          end
          unless m[3].nil?
            type = ["null", type]
          end

          return type
        end

        val
      end

      def load(val)
        if val.is_a?(Array)
          @inner.load(val.map do |item|
            item.is_a?(String) ? resolve(item) : item
          end)
        elsif val.is_a?(String)
          @inner.load(resolve(val))
        else
          @inner.load(val)
        end
      end

      OUTER_TYPE_LOADER = TypedDSLLoader.new(
          UnionLoader.new(
              RecordType::RecordTypeLoader.new,
              ArrayType::InnerLoader.new,
              EnumType::InnerLoader.new,
              EnumLoader::PRIMITIVE_TYPE,
              EnumLoader::NOMINAL_TYPE,
          )
      )

      WITH_UNIONS_TYPE_LOADER = TypedDSLLoader.new(
          UnionLoader.new(
              OUTER_TYPE_LOADER,
              ArrayLoader.new(OUTER_TYPE_LOADER),
          )
      )
    end

    def self.load_item(val, field_type)
      if val.is_a?(Hash)
        if val.include?("$import")
          raise "$import expressions are not yet supported"
        elsif val.include?("$include")
          raise "$include expressions are not yet supported"
        end
      end

      return field_type.load(val)
    end

    class InputParameter < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          label: PrimitiveLoader::STRING.optional,
          secondaryFiles: NeverLoader::UNSUPPORTED,
          streamable: NeverLoader::UNSUPPORTED,
          loadContents: NeverLoader::UNSUPPORTED,
          loadListing: NeverLoader::UNSUPPORTED,
          valueFrom: NeverLoader::UNSUPPORTED,
          doc: PrimitiveLoader::STRING.optional,

          type: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
          default: AnyLoader::ANY.optional,
          format: PrimitiveLoader::STRING.optional,
      }

      def default
        default = @attributes['default']
        return nil unless default
        RecordType::Field.type_loader(@attributes['type'])&.load(default)
      end
    end

    class OutputParameter < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          label: PrimitiveLoader::STRING.optional,
          secondaryFiles: NeverLoader::UNSUPPORTED,
          streamable: NeverLoader::UNSUPPORTED,
          doc: PrimitiveLoader::STRING.optional,

          outputBinding: NeverLoader::UNSUPPORTED,
          type: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
          format: PrimitiveLoader::STRING.optional,
      }
    end

    class WorkflowOutputParameter < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING,
          label: PrimitiveLoader::STRING.optional,
          secondaryFiles: NeverLoader::UNSUPPORTED,
          streamable: NeverLoader::UNSUPPORTED,
          linkMerge: NeverLoader::UNSUPPORTED,
          pickValue: NeverLoader::UNSUPPORTED,
          doc: PrimitiveLoader::STRING.optional,

          outputSource: SourceLoader.new,
          type: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
          format: PrimitiveLoader::STRING.optional,
      }

      def outputSource
        @attributes['outputSource']
      end

      def id
        @attributes['id']
      end

      def type
        @attributes['type']
      end

      def format
        @attributes['format']
      end
    end

    class WorkflowInputParameter < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING,
          label: PrimitiveLoader::STRING.optional,
          secondaryFiles: NeverLoader::UNSUPPORTED,
          streamable: NeverLoader::UNSUPPORTED,
          loadContents: NeverLoader::UNSUPPORTED,
          loadListing: NeverLoader::UNSUPPORTED,
          doc: PrimitiveLoader::STRING.optional,
          inputBinding: NeverLoader::UNSUPPORTED,

          default: AnyLoader::ANY,
          type: TypedDSLLoader::WITH_UNIONS_TYPE_LOADER,
          format: PrimitiveLoader::STRING.optional,
      }

      def id
        @attributes['id']
      end

      def default
        @attributes['default']
      end

      def type
        @attributes['type']
      end
    end

    class StepOutput < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING,
      }

      def id
        @attributes['id']
      end
    end


    class StepInput < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          source: SourceLoader.new.optional,
          label: PrimitiveLoader::STRING.optional,
          linkMerge: NeverLoader::UNSUPPORTED,
          pickValue: NeverLoader::UNSUPPORTED,
          loadContents: NeverLoader::UNSUPPORTED,
          loadListing: NeverLoader::UNSUPPORTED,
          valueFrom: NeverLoader::UNSUPPORTED,
          default: AnyLoader::ANY.optional,
      }

      def id
        @attributes['id']
      end

      def source
        @attributes['source']
      end
    end

    class Operation < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          label: PrimitiveLoader::STRING.optional,
          doc: PrimitiveLoader::STRING.optional,
          requirements: NeverLoader::UNSUPPORTED,
          hints: NeverLoader::UNSUPPORTED,
          cwlVersion: EnumLoader.new("v1.0", "v1.1", "v1.2").optional,
          intent: NeverLoader::UNSUPPORTED,
          class: EnumLoader.new("Operation"),
          inputs: InputParameter.loader.as_mapped_array('id', 'type'),
          outputs: OutputParameter.loader.as_mapped_array('id', 'type'),
      }

      def id
        @attributes['id']
      end
    end

    class Step < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          label: PrimitiveLoader::STRING.optional,
          doc: PrimitiveLoader::STRING.optional,
          in: StepInput.loader.as_mapped_array('id', 'source'),
          out: StepOutput.loader.or(PrimitiveLoader::STRING.map { |id| StepOutput.loader.load({'id' => id}) }).as_array,
          requirements: NeverLoader::UNSUPPORTED,
          hints: NeverLoader::UNSUPPORTED,
          run: PrimitiveLoader::STRING.map { |id| Operation.loader.load({'id' => id, 'class' => 'Operation', 'inputs' => [], 'outputs' => []}) }.or(Operation.loader),
          when: NeverLoader::UNSUPPORTED,
          scatter: NeverLoader::UNSUPPORTED,
          scatterMethod: NeverLoader::UNSUPPORTED,
      }

      def id
        @attributes['id']
      end

      def in
        @attributes['in']
      end

      def out
        @attributes['out']
      end
    end

    class Workflow < Cwl
      FIELD_LOADERS = {
          id: PrimitiveLoader::STRING.optional,
          label: PrimitiveLoader::STRING.optional,
          doc: PrimitiveLoader::STRING.optional,
          requirements: NeverLoader::UNSUPPORTED,
          hints: NeverLoader::UNSUPPORTED,
          intent: NeverLoader::UNSUPPORTED,
          class: EnumLoader.new("Workflow"),
          cwlVersion: EnumLoader.new("v1.0", "v1.1", "v1.2"),
          inputs: WorkflowInputParameter.loader.as_mapped_array('id', 'type'),
          outputs: WorkflowOutputParameter.loader.as_mapped_array('id', 'type'),
          steps: Step.loader.as_mapped_array('id', 'source')
      }

      def inputs
        @attributes['inputs']
      end

      def outputs
        @attributes['outputs']
      end

      def steps
        @attributes['steps']
      end
    end
  end
end