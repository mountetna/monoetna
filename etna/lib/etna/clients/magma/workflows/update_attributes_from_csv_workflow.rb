require 'csv'
require 'ostruct'
require_relative './crud_workflow'

module Etna
  module Clients
    class Magma
      class UpdateAttributesFromCsvWorkflowBase < Struct.new(:magma_crud, :project_name, :filepath, :model_name, keyword_init: true)
        def initialize(opts)
          super(**{}.update(opts))
        end

        def magma_client
          magma_crud.magma_client
        end

        def parse_input_file
          raise "Must be implemented in a subclass"
        end

        def model_exists?(model_name)
          models.model(model_name)
        end

        def consolidate_attributes_to_hash(row)
          raise "Must be implemented in a subclass"
        end

        def find_attribute(model_name, attribute_name)
          models.model(model_name).template.attributes.attribute(attribute_name)
        end

        def each_revision
          parse_input_file do |model_name, record_name, revision|
            yield [model_name, record_name, revision]
          end
        end

        def update_attributes
          magma_crud.update_records do |update_request|
            each_revision do |model_name, record_name, revision|
              update_request.update_revision(model_name, record_name, revision)
            end
          end
        end

        def models
          @models ||= begin
            magma_client.retrieve(RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
          end
        end
      end

      class RowBase
        def stripped_value(attribute_value)
          attribute_value ? attribute_value.strip : attribute_value
        end

        def nil_or_empty?(value)
          value.nil? || value.empty?
        end
      end

      class UpdateAttributesFromCsvWorkflowMultiModel < UpdateAttributesFromCsvWorkflowBase
        def parse_input_file
          CSV.parse(File.read(filepath)).map do |csv_row|
            row = Row.new(csv_row, self)

            raise "Invalid model \"#{row.model_name}\" for project #{project_name}." unless model_exists?(row.model_name)

            yield [row.model_name, row.record_name, row.to_h]
          end
        end

        class Row < RowBase
          attr_reader :model_name, :record_name
          def initialize(raw, workflow)
            # Assumes rows are in pairs, where
            #   [0] = model_name
            #   [1] = record_name / identifier
            #   [2], [4], etc. = attribute_name
            #   [3], [5], etc. = attribute value
            # So not every row needs the same number of columns
            @raw = raw
            @workflow = workflow

            raise "Invalid revision row #{@raw}. Must include at least 4 column values (model,record_name,attribute_name,attribute_value)." if @raw.length < 4
            raise "Invalid revision row #{@raw}. Must have an even number of columns." if @raw.length.odd?

            @model_name = raw[0]

            raise "Invalid model name: \"#{@model_name}\"." if nil_or_empty?(@model_name)

            @model_name.strip!

            @record_name = raw[1]

            raise "Invalid record name: \"#{@record_name}\"." if nil_or_empty?(@record_name)

            @record_name.strip!
          end

          def to_h
            # Take attribute index values (even) and put them into a hash.
            # The assigned value will be the subsequent odd index.
            # {attribute_name: attribute_value}
            {}.tap do |attributes|
              (2..(@raw.length - 1)).to_a.each do |index|
                if index.even?
                  attribute_name = @raw[index]

                  raise "Invalid attribute name: \"#{attribute_name}\"." if nil_or_empty?(attribute_name)
                  attribute_name.strip!

                  raise "Invalid attribute #{attribute_name} for model #{model_name}." unless attribute = @workflow.find_attribute(model_name, attribute_name)

                  attributes[attribute_name] = stripped_value(@raw[index + 1])
                end
              end
            end
          end
        end
      end

      class UpdateAttributesFromCsvWorkflowSingleModel < UpdateAttributesFromCsvWorkflowBase
        def initialize(opts)
          super(**{}.update(opts))
          raise "Single Model invokation must include keyword :model_name." if !opts[:model_name]
          raise "Invalid model #{model_name} for project #{project_name}." unless model_exists?(model_name)
        end

        def parse_input_file
          CSV.parse(File.read(filepath), headers: true).map do |csv_row|
            row = Row.new(csv_row, model_name, self)

            yield [model_name, row.record_name, row.to_h]
          end
        end

        class Row < RowBase
          attr_reader :record_name
          def initialize(raw, model_name, workflow)
            # Assumes CSV includes a column header to identify the attribute_name
            # Assumes index 0 is the record_name
            @raw = raw
            @model_name = model_name
            @workflow = workflow

            @record_name = @raw[0]
            raise "Invalid record name: \"#{record_name}\"." if nil_or_empty?(record_name)

            @record_name.strip!
          end

          def to_h
            # Row can be converted to a hash, where keys are attribute_names and the
            #   values come from the CSV
            # {attribute_name: attribute_value}
            {}.tap do |attributes|
              row_hash = @raw.to_h
              row_keys = row_hash.keys
              row_keys[1..row_keys.length - 1].each do |attribute_name|

                raise "Invalid attribute name: \"#{attribute_name}\"." if nil_or_empty?(attribute_name)

                attribute_name_clean = attribute_name.strip
                raise "Invalid attribute \"#{attribute_name_clean}\" for model #{@model_name}." unless attribute = @workflow.find_attribute(@model_name, attribute_name_clean)

                attributes[attribute_name_clean] = stripped_value(@raw[attribute_name])
              end
            end
          end
        end
      end
    end
  end
end

