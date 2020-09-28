require 'csv'
require 'ostruct'
require_relative './crud_workflow'

module Etna
  module Clients
    class Magma
      class UpdateAttributesFromCsvWorkflow < Struct.new(:magma_crud, :project_name, :filepath, keyword_init: true)
        def initialize(opts)
          super(**{}.update(opts))
        end

        def magma_client
          magma_crud.magma_client
        end

        def parse_input_file
          CSV.parse(File.read(filepath)).map do |row|
            # Assumes rows are in pairs, where
            #   [0] = model_name
            #   [1] = record_name / identifier
            #   [2], [4], etc. = attribute_name
            #   [3], [5], etc. = attribute value
            # So not every row needs the same number of columns
            raise "Invalid model #{row[0]} for project #{project_name}." unless model_exists?(row[0])

            yield [row[0], row[1], consolidate_attributes_to_hash(row)]
          end
        end

        def model_exists?(model_name)
          models.model(model_name)
        end

        def consolidate_attributes_to_hash(row)
          # Take attribute index values (even) and put them into a hash.
          # The assigned value will be the subsequent odd index. Make sure
          #   to cast it to the right format, since CSV may treat everything
          #   as a string.
          # {attribute_name: attribute_value}
          {}.tap do |attributes|
            (2..(row.length - 1)).to_a.each do |index|
              if index % 2 == 0
                attribute_name = row[index]

                next if attribute_name.empty?
                raise "Invalid attribute #{attribute_name} for model #{row[0]}." unless attribute = find_attribute(row[0], attribute_name)

                attribute_value = format_value(attribute, row[index + 1])

                attributes[attribute_name] = attribute_value
              end
            end
          end
        end

        def find_attribute(model_name, attribute_name)
          models.model(model_name).template.attributes.attribute(attribute_name)
        end

        def format_value(attribute, attribute_value)
          return nil unless attribute_value

          # Make sure to cast the value from the CSV
          # More complex types not currently supported
          case attribute.attribute_type
          when Etna::Clients::Magma::AttributeType::BOOLEAN
            attribute_value.to_s.downcase == 'true' || attribute_value.to_s == '1'
          when Etna::Clients::Magma::AttributeType::INTEGER
            attribute_value.to_i
          when Etna::Clients::Magma::AttributeType::FLOAT
            attribute_value.to_f
          else
            attribute_value
          end
        end

        def each_revision
          parse_input_file do |model_name, record_name, revision|
            yield [model_name, record_name, revision]
          end
        end

        def update_attributes
          # Want to send JSON to the /update route, so
          #   use magma_client directly.
          each_revision do |model_name, id, revision|
            update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project_name)
            update_request.update_revision(model_name, id, revision)
            magma_client.update(update_request, true)
          end
        end

        def models
          @models ||= begin
            magma_client.retrieve(RetrievalRequest.new(project_name: self.project_name, model_name: 'all')).models
          end
        end
      end
    end
  end
end

