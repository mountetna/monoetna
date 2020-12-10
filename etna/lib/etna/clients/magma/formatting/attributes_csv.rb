require 'ostruct'
require_relative '../../../csvs'

module Etna
  module Clients
    class Magma
      module AttributesCsv
        class Exporter < Etna::CsvExporter
          attr_reader :attributes_to_export_by_model_name

          def initialize(attributes_to_export_by_model_name = {})
            @attributes_to_export_by_model_name = attributes_to_export_by_model_name
            max_attributes = attributes_to_export_by_model_name.map { |k, attributes| attributes.length }.max || 0

            super([
                :model_name,
                :identifier,
                :table_attribute,
            ] + max_attributes.times.map { |i| "attribute_#{i}" })
          end

          def write_attributes(models, documents, output_io: nil, filename: nil)
            with_row_writeable(filename: filename, output_io: output_io) do |row_writeable|
              write_attribute_rows(models, documents, row_writeable)
            end
          end

          def write_attribute_rows(models, documents, row_writeable)
            @attributes_to_export_by_model_name.each do |model_name, attributes|
              self.write_model_attribute_rows(models, model_name, documents, row_writeable)
            end
          end

          def header_row(model, **base)
            row_attributes, _ = self.row_and_table_attributes(model)
            base.update(row_attributes.each_with_index.map { |attr_name, i| ["attribute_#{i}", attr_name] }.to_h)
          end

          def write_model_attribute_rows(models, model_name, documents, row_writeable)
            return unless (model = models.model(model_name))
            row_writeable.write
            row_writeable << header_row(model, model_name: model_name)
            row_attributes, table_attributes = self.row_and_table_attributes(model_name)

            model_docs = documents[model_name] || {}
            model_docs.each do |ident, doc|
              row_writeable << ({identity: ident}.update(row_attributes.map do |attr|
                val = doc[attr.name]
                self.class.serialize_attr(attr, val)
              end))

              table_attributes.each do |table_attr|
                row_writeable << header_row(models.model(table_attr.link_model_name), table_attribute: table_attr.name)

                doc[table_attr.name].each do |table_ident|
                  table_row_doc = documents[table_attr.link_model_name][table_ident]
                  table_row_attributes, _ = row_and_table_attributes(table_attr.link_model_name)

                  row_writeable.write(table_row_attributes.map do |table_row_attr|
                    val = table_row_doc[table_row_attr.name]
                    self.class.serialize_attr(table_row_attr, val)
                  end)
                end
              end
            end
          end

          def row_and_table_attributes(model)
            model_export_attribute_names = @attributes_to_export_by_model_name[model.name] || []
            row_attributes = model_export_attribute_names.map { |attr_name| model.template.attributes.attribute(attr_name) }
            table_attributes = row_attributes.select { |attr| attr.attribute_type == Etna::Clients::Magma::AttributeType::TABLE }
            row_attributes = row_attributes.select { |attr| !table_attributes.include?(attr) }

            [row_attributes, table_attributes]
          end

          def self.serialize_attr(attr, val)
            case attr.attribute_type
            when AttributeType::LINK, AttributeType::IDENTIFIER, AttributeType::PARENT, AttributeType::STRING
              val || ''
            when AttributeType::TABLE
              raise ImportError.new(format_attr_err_message(attr, "cannot be set in a single row column, use the table_attribute column instead."))
            when AttributeType::COLLECTION, AttributeType::CHILD
              raise TypeError.new(format_attr_err_message(attr, "cannot be set in a single row column", recip: recip))
            when AttributeType::MATRIX
              # TODO: Make this an external file loading?
              JSON.stringify(val)
            when AttributeType::FILE_COLLECTION
              if val.empty?
                ""
              else
                JSON.stringify(val.map { |rev| serialize_file(rev) })
              end
            when AttributeType::INTEGER, AttributeType::FLOAT
              val.to_s
            when AttributeType::BOOLEAN
              val ? 'true' : 'false'
            when AttributeType::DATE_TIME
              DateTime.parse(val).iso8601
            when AttributeType::MATCH
              JSON.stringify(val)
            when AttributeType::FILE
              serialize_file(val)
            else
              raise TypeError.new("#{attr.attribute_type} is not supported in this version.")
            end
          end

          def self.serialize_file(rev)
            result = rev['path']
            if (original_filename = rev['original_filename'])
              result += ",#{original_filename}"
            end
            result
          end
        end

        class Importer < Etna::CsvImporter
          def initialize(models)
            @models = models
            super()
          end

          def process_row(row_processor, update_request, models)
            variable_keys = row_processor.row.keys.select { |k| k.start_with?('attribute_') }

            row_processor.process(:model_name) do |model_name|
              model = models.model(model_name)
              raise ImportError.new("model #{model_name} does not exist in the target project.") if model.nil?
              [model, variable_keys.map { |k| row_processor.row[k] }]
            end

            row_processor.process(:identifier, :model_name) do |identifier, model_and_columns|
              model, columns = model_and_columns
              values = variable_keys.zip(row).map { |k, v| [columns[k], v] }.to_h
              update_request.update_revision(model.name, identifier, prepare_update(models, model, values.update({model.template.identifier => identifier})))
              identifier
            end

            unless row_processor.process(:table_attribute, :model_name, :identifier) do |table_attribute_name, model_and_model_columns, identifier|
              model, _ = model_and_model_columns
              columns = variable_keys.map { |k| row_processor.row[k] }

              unless (attribute = model.template.attributes.attribute(table_attribute_name))
                raise ImportError.new("attribute #{table_attribute_name} does not exist on the #{model.name} model.")
              end

              unless attribute.attribute_type == Etna::Clients::Magma::AttributeType::TABLE
                raise ImportError.new("attribute #{table_attribute_name} on the #{model.name} model is not a table")
              end

              unless (table_model = models.model(attribute.link_model_name))
                raise ImportError.new("Badly configured source models; cannot find model named #{attribute.link_model_name}")
              end

              ->(row) do
                values = variable_keys.zip(row).map { |k, v| [columns[k], v] }.to_h
                update_request.append_table(model.name, identifier, table_attribute_name, prepare_update(models, table_model, values))
              end
            end
              if (table_appender = row_processor.context[:table_attribute])
                variable_row = variable_keys.map { |k| row_processor.row[k] }
                table_appender.call(variable_row)
              end
            end
          end

          def prepare_update_request(models, model_name: nil, filename: nil, io: nil, &validation_err_block)
            Etna::Clients::Magma::UpdateRequest.new.tap do |result|
              context = {}
              each_csv_row(filename: filename, input_io: io) do |row, lineno|
                unless model_name.nil?
                  row[:model_name] = model_name
                end

                process_row(Etna::CsvImporter::NestedRowProcessor.new(row, lineno, context), result, models)
              end
            end
          rescue ImportError => e
            validation_err_block.call(e.message)
          end

          def prepare_update(models, model, update_hash)
            update_hash.keys.each do |attr_name|
              unless (attr = model.template.attributes.attribute(attr_name))
                raise ImportError.new(format_attr_err_message(attr, "does not exist!"))
              end

              recip = models.find_reciprocal(model: model, link_attribute_name: attr_name)
              string_value = update_hash[attr_name].strip

              update_hash[attr_name] = self.class.deserialize(attr, string_value, recip)
            end
          end

          METIS_PATH_REGEX = Regexp.new('^metis://' +
              # project_name
              '([^/]*?)/' +
              # bucket_name
              '([^/]*?)/' +
              # folder path + filename
              '(.*)$')

          def self.deserialize(attr, string_value, recip)
            case attr.attribute_type
            when AttributeType::LINK, AttributeType::IDENTIFIER, AttributeType::PARENT, AttributeType::STRING
              string_value
            when AttributeType::TABLE
              raise ImportError.new(format_attr_err_message(attr, "cannot be set in a single row column, use the table_attribute column instead."))
            when AttributeType::COLLECTION, AttributeType::CHILD
              raise TypeError.new(format_attr_err_message(attr, "cannot be set in a single row column", recip: recip))
            when AttributeType::MATRIX
              begin
                # TODO: Make this an external file loading?
                JSON.parse(string_value)
              rescue
                raise ImportError.new(format_attr_err_message(attr, "must be a json encoded array of numbers"))
              end
            when AttributeType::FILE_COLLECTION
              if string_value.empty?
                []
              else
                begin
                  JSON.parse(update_hash).map { |s| self.prepare_file_value(s) }
                rescue JSONError
                  raise ImportError.new(format_attr_err_message(attr, 'must be a valid json encoded array of strings.'))
                end
              end
            when AttributeType::INTEGER
              update_hash[attr_name] = to_numeric_or_error(string_value, [Integer])
            when AttributeType::BOOLEAN
              if ['true', 'y', 't', 'yes', 'false', 'n', 'f', 'no'].include?(string_value.downcase)
                Etna::CsvImporter::COLUMN_AS_BOOLEAN.call(string_value)
              else
                raise ImportError.new(format_attr_err_message(attr, "must be one of true, false, yes, no"))
              end
            when AttributeType::DATE_TIME
              begin
                DateTime.parse(string_value).iso8601
              rescue Date::Error
                raise ImportError.new(format_attr_err_message(attr, "is not supported"))
              end
            when AttributeType::MATCH
              JSON.parse(string_value)
            when AttributeType::FLOAT
              to_numeric_or_error(string_value, Float)
            when AttributeType::FILE
              self.prepare_file_value(string_value)
            else
              raise ImportError.new(format_attr_err_message(attr, "is not supported"))
            end
          end
        end

        def self.prepare_file_value(string_value)
          if string_value.empty?
            {path: '::blank'}
          else
            original_filename = nil
            comma_idx = string_value.index(',')
            if comma_idx
              original_filename = string_value.slice((comma_idx + 1)..-1)
              string_value = string_value.slice(0..comma_idx)
            end

            if string_value =~ METIS_PATH_REGEX
              if original_filename.nil?
                original_filename = File.basename(string_value)
              end
              {path: string_value, original_filename: original_filename}
            else
              raise ArgumentError.new("must be a valid metis path (metis://project/bucket/path)")
            end
          end
        end

        def self.format_attr_err_message(attr, msg, recip: nil)
          msg = "#{attr.attribute_type} attribute #{attr.name} #{msg}"
          if recip
            msg += ", try using the reciprocal attribute #{recip.name} on #{attr.link_model_name}"
          end

          msg
        end

        def self.to_numeric_or_error(string_value, numeric_types = [Integer, Float])
          numeric_types.each do |type|
            begin
              return type.call(string_value)
            rescue ArgumentError => e
              next
            end
          end

          raise ArgumentError.new("must be a valid #{numeric_types.map(&:name).join(' or ')}")
        end
      end
    end
  end
end
