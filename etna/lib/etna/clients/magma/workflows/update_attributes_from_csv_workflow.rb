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
            raise "Invalid revision row #{row}. Must include at least 4 column values (model,record_name,attribute_name,attribute_value)." if row.length < 4
            raise "Invalid revision row #{row}. Must have an even number of columns." if row.length % 2 == 1

            model_name = row[0]

            raise "Invalid model #{model_name} for project #{project_name}." unless model_exists?(model_name)

            yield [model_name, row[1], consolidate_attributes_to_hash(row)]
          end
        end

        def model_exists?(model_name)
          models.model(model_name)
        end

        def consolidate_attributes_to_hash(row)
          # Take attribute index values (even) and put them into a hash.
          # The assigned value will be the subsequent odd index.
          # {attribute_name: attribute_value}
          {}.tap do |attributes|
            (2..(row.length - 1)).to_a.each do |index|
              if index % 2 == 0
                attribute_name = row[index]

                next if attribute_name.empty?
                raise "Invalid attribute #{attribute_name} for model #{row[0]}." unless attribute = find_attribute(row[0], attribute_name)

                attributes[attribute_name] = row[index + 1]
              end
            end
          end
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
    end
  end
end

