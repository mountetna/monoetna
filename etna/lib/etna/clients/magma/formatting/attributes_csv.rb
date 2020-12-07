require 'ostruct'
require_relative '../../../csvs'

module Etna
  module Clients
    class Magma
      module AttributesCsv
        class Importer < Etna::CsvImporter
          def initialize(models)
            @models = models
            super(&method(:format_row))
          end

          def format_row(row)
            replace_row_column(row, :attribute_type) { |s| AttributeType.new(s) }
            replace_row_column(row, :desc) { row.delete(:description) }
            replace_row_column(row, :restricted, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :read_only, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :options) { |s| {"type" => "Array", "value" => s.split(',').map(&:strip)} }
            replace_row_column(row, :match) { |s| {"type" => "Regexp", "value" => Regexp.new(s).source} }
            replace_row_column(row, :validation) { row[:options] || row[:match] }
            replace_row_column(row, :hidden, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :unique, &COLUMN_AS_BOOLEAN)
          end

          def process_row(row_processor = NestedRowProcessor.new, update_request = Etna::Clients::Magma::UpdateRequest.new, models = Etna::Clients::Magma::Models.new)
            row_processor.process(:model_name) do |model_name|
              model = models.model(model_name)
              raise ImportError.new("model #{model_name} does not exist in the target project.") if model.nil?
              model
            end

            row_processor.process(:identifier, :model_name)
            row_processor.row.keys.select { |k| k != :model_name && k != :identifier }.each do |col_name|
              row_processor.watch(col_name, :model_name, :identifier) do |prev_col_val, col_val, model, identifier|
                if prev_col_val.nil?
                  next col_val
                end

                attr_name = prev_col_val
                unless (attr = model.template.attributes.attribute(attr_name))
                  raise ImportError.new("Attribute #{attr_name} does not exist on model #{model.name}")
                end

                recip = models.find_reciprocal(model: model, link_attribute_name: attr)
                attr_value = Etna::Clients::Magma.cast(col_val, attr, recip)
                next prev_col_val
              end
            end
          end

          def prepare_update_request(models, filename: nil, io: nil, &validation_err_block)
            update_request = Etna::Clients::Magma::UpdateRequest.new

            each_csv_row(filename: filename, input_io: io) do |row, lineno|
              p = NestedRowProcessor.new(row, lineno, context)
              process_row(p, update_request, models)
            end

            update_request
          rescue ImportError => e
            validation_err_block.call(e.message)
          end
        end
      end
    end
  end
end
