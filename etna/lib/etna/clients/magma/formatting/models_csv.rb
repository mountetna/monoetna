require 'ostruct'
require_relative '../../../csvs'

module Etna
  module Clients
    class Magma
      module ModelsCsv
        module Prettify
          def prettify(name)
            name.split('_').map(&:capitalize).join(' ')
          end
        end

        class Exporter < Etna::CsvExporter
          include Prettify

          def initialize
            super([
                :comments,
                :model_name, :identifier, :parent_model_name, :parent_link_type,
                :attribute_name,
                :new_attribute_name,
                :attribute_type,
                :link_model_name,
                :description,
                :display_name,
                :format_hint,
                :restricted,
                :read_only,
                :options,
                :match,
                :attribute_group,
                :hidden,
                :unique,
                :matrix_constant,
            ])
          end

          def write_models(models, model_keys = models.model_keys.sort, output_io: nil, filename: nil)
            with_row_writeable(filename: filename, output_io: output_io) do |row_writeable|
              write_model_rows(models, model_keys, row_writeable)
            end
          end

          def write_model_rows(models, model_keys, row_writeable)
            ensure_parents(models, model_keys, row_writeable)

            matrix_constants = {}
            model_keys.each { |model_name| each_model_row(models, model_name, matrix_constants, row_writeable) }

            matrix_constants.each do |digest, options|
              row_writeable.write
              row_writeable.write(matrix_constant: COPY_OPTIONS_SENTINEL + digest)
              options.each do |option|
                row_writeable.write(matrix_constant: option)
              end
            end
          end

          def ensure_parents(models, model_keys, row_writeable)
            q = model_keys.dup
            seen = Set.new

            until q.empty?
              model_key = q.shift
              next if model_key.nil?
              next if seen.include?(model_key)
              seen.add(model_key)

              # For models that are only part of the trunk, but not of the tree of model_keys,
              # we still need their basic information (identifier / parent) for validation and for
              # potentially creating the required tree dependencies to connect it to a remote tree.
              unless model_keys.include?(model_key)
                each_model_trunk_row(models, model_key, row_writeable)
              end

              q.push(*models.model(model_key).template.all_linked_model_names)
            end
          end

          def each_model_trunk_row(models, model_name, row_writeable)
            return unless (model = models.model(model_name))

            # Empty link for better visual separation
            row_writeable.write
            row_writeable.write(model_name: model_name)

            if !model.template.parent.nil?
              parent_model = models.model(model.template.parent)
              reciprocal = models.find_reciprocal(model: model, link_model: parent_model)

              row_writeable.write(
                  identifier: model.template.identifier,
                  parent_model_name: model.template.parent,
                  parent_link_type: reciprocal.attribute_type.to_s
              )
            else
              row_writeable.write(
                  identifier: model.template.identifier,
              )
            end
          end

          def each_model_row(models, model_name, matrix_constants, row_writeable)
            return unless (model = models.model(model_name))

            each_model_trunk_row(models, model_name, row_writeable)
            model.template.attributes.all.each do |attribute|
              each_attribute_row(models, model, attribute, matrix_constants, row_writeable)
            end
          end

          def each_attribute_row(models, model, attribute, matrix_constants, row_writeable)
            if attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
              # Identifiers for models whose parent link type ends up being a table are non configurable, so we don't
              # want to include them in the CSV.
              if models.find_reciprocal(model: model, link_attribute_name: model.template.parent)&.attribute_type == Etna::Clients::Magma::AttributeType::TABLE
                return
              end
            else
              return unless Etna::Clients::Magma::AttributeValidator.valid_add_row_attribute_types.include?(attribute.attribute_type)
            end

            options = attribute.options&.join(', ')
            if attribute.attribute_type == Etna::Clients::Magma::AttributeType::MATRIX
              # Matrix attribute validations are massive, and functional.  For now, we don't support showing and editing
              # them inline with this spreadsheet.  In the future, I think we should possibly introduce the concept of
              # CONSTANTS or Matrix Types that are managed separately.
              options = options || ''
              digest = Digest::MD5.hexdigest(options)
              matrix_constants[digest] ||= options.split(',').map(&:strip)

              options = COPY_OPTIONS_SENTINEL + digest
            end

            row_writeable.write(
                attribute_name: attribute.name,
                attribute_type: attribute.attribute_type,
                link_model_name: attribute.link_model_name,
                reciprocal_link_type: models.find_reciprocal(model: model, attribute: attribute)&.attribute_type,
                description: attribute.desc,
                display_name: attribute.display_name,
                match: attribute.match,
                format_hint: attribute.format_hint,
                restricted: attribute.restricted,
                read_only: attribute.read_only,
                options: options,
                attribute_group: attribute.attribute_group,
                hidden: attribute.hidden,
                unique: attribute.unique,
            )
          end
        end

        class ModelsChangeset < Struct.new(:models, :renames, :matrix_constants, keyword_init: true)
          def initialize(*args)
            super

            self.models ||= Models.new
            self.renames ||= {}
            self.matrix_constants ||= {}
          end

          def build_renames(model_name)
            renames[model_name] ||= {}
          end
        end

        COPY_OPTIONS_SENTINEL = '$copy_from_original$'

        class Importer < Etna::CsvImporter
          include Prettify

          # Columns of the row, _after format_row is called_, that should be applied to an attribute.
          # This should line up with the attribute names _on the model itself_.
          ATTRIBUTE_ROW_ENTRIES = [
              :attribute_type,
              :link_model_name, :description,
              :display_name, :format_hint,
              :restricted, :read_only,
              :validation, :attribute_group,
              :hidden, :unique,
          ]

          def initialize
            super(&method(:format_row))
          end

          def format_row(row)
            replace_row_column(row, :attribute_type) { |s| AttributeType.new(s) }
            # replace_row_column(row, :desc) { row.delete(:description) }
            replace_row_column(row, :restricted, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :read_only, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :options) { |s| {"type" => "Array", "value" => s.split(',').map(&:strip)} }
            replace_row_column(row, :match) { |s| {"type" => "Regexp", "value" => Regexp.new(s).source} }
            replace_row_column(row, :validation) { row[:options] || row[:match] }
            replace_row_column(row, :hidden, &COLUMN_AS_BOOLEAN)
            replace_row_column(row, :unique, &COLUMN_AS_BOOLEAN)
          end

          def process_row(row_processor, changeset)
            models = changeset.models

            process_matrix_constants(changeset, row_processor)
            process_model_config(models, row_processor)
            process_parent_config(models, row_processor)
            process_attribute_properties(changeset, row_processor)
          end

          def prepare_changeset(filename: nil, input_io: nil, &validation_err_block)
            ModelsChangeset.new.tap do |changeset|
              context = {}
              each_csv_row(filename: filename, input_io: input_io) do |row, lineno|
                p = NestedRowProcessor.new(row, lineno, context)
                process_row(p, changeset)
              end

              # After configuring attributes, set defaults for certain attribute properties.
              changeset.models.all.each do |model|
                model.template.attributes.all.each do |attribute|
                  attribute.set_field_defaults!
                end
              end
            end
          rescue ImportError => e
            validation_err_block.call(e.message)
          end


          private

          def process_model_config(models, row_processor)
            row_processor.process(:model_name) do |model_name|
              models.build_model(model_name).build_template.tap { |t| t.name = model_name }
            end
            row_processor.process(:identifier, :model_name) do |identifier, template|
              template.identifier = identifier
            end
          end

          def process_parent_config(models, row_processor)
            row_processor.process(:parent_model_name, :model_name) do |parent_model_name, template|
              if template.parent && !template.parent.empty? && template.parent != parent_model_name
                raise ImportError.new("Model #{template.name} was provided multiple parent_model_names: #{template.parent} and #{parent_model_name}")
              end

              template.parent = parent_model_name

              template.build_attributes.build_attribute(template.parent).tap do |parent_att|
                parent_att.name = parent_att.attribute_name = parent_model_name
                parent_att.attribute_type = Etna::Clients::Magma::AttributeType::PARENT
                parent_att.link_model_name = parent_model_name
                parent_att.description = prettify(parent_model_name)
                parent_att.display_name = prettify(parent_model_name)
              end
            end

            row_processor.process(:parent_link_type, :model_name, :parent_model_name) do |parent_link_type, template|
              reciprocal = models.find_reciprocal(model_name: template.name, link_attribute_name: template.parent)
              if reciprocal && reciprocal.attribute_type.to_s != parent_link_type
                raise ImportError.new("Model #{template.name} was provided multiple parent_link_types: #{reciprocal.attribute_type} and #{parent_link_type}")
              end

              if reciprocal && reciprocal.attribute_name != template.name
                raise ImportError.new("Model #{template.name} is linked to #{template.parent}, but the reciprocal link is misnamed as '#{reciprocal.attribute_name}'.")
              end

              models.build_model(template.parent).tap do |parent_model|
                parent_model_attribute_by_model_name = parent_model.build_template.build_attributes.attribute(template.name)
                if parent_model_attribute_by_model_name && !reciprocal
                  raise ImportError.new("Model #{template.parent} is linked as a parent to #{template.name}, but it already has an attribute named #{template.name} #{parent_model_attribute_by_model_name.raw}.")
                end

                parent_model.build_template.build_attributes.build_attribute(template.name).tap do |attr|
                  attr.attribute_name = attr.name = template.name
                  attr.attribute_type = parent_link_type
                  attr.link_model_name = template.name
                  attr.description = prettify(template.name)
                  attr.display_name = prettify(template.name)
                end
              end
            end
          end

          def process_matrix_constants(changeset, row_processor)
            row_processor.process(:matrix_constant) do |matrix_constant|
              # Lines that start with the sentinel are a distinct 'group'
              if matrix_constant.start_with?(COPY_OPTIONS_SENTINEL)
                matrix_constant = matrix_constant.slice((COPY_OPTIONS_SENTINEL.length)..-1)
                changeset.matrix_constants[matrix_constant] = []
              else
                # Until the sentinel is seen again, all further entries belong to the set.
                changeset.matrix_constants[changeset.matrix_constants.keys.last] << matrix_constant
              end

              matrix_constant
            end
          end

          def process_attribute_properties(changeset, row_processor)
            models = changeset.models

            row_processor.process(:attribute_name, :model_name) do |attribute_name, template|
              attributes = template.build_attributes

              existing_attribute = attributes.attribute(attribute_name)
              if existing_attribute
                if existing_attribute.attribute_type != Etna::Clients::Magma::AttributeType::COLLECTION
                  raise ImportError.new("Attribute #{attribute_name} of model #{template.name} has duplicate definitions!")
                end

                # Clear its definition; implicit attributes are overwritten by explicit definitions if they exist.
                attributes.raw.delete(attribute_name)
              end

              attributes.build_attribute(attribute_name).tap { |att| att.name = att.attribute_name = attribute_name }
            end

            row_processor.process(:new_attribute_name, :attribute_name, :model_name) do |new_attribute_name, att, template|
              renames = changeset.build_renames(template.name)
              if renames.include?(att.name) && renames[att.name] != new_attribute_name
                raise ImportError.new("Found multiple new_attribute_name values for #{template.name}'s #{att.name}': #{new_attribute_name} or #{renames[att.name]}?")
              end

              renames[att.attribute_name] = new_attribute_name
            end

            ATTRIBUTE_ROW_ENTRIES.each do |prop_name|
              row_processor.process(prop_name, :model_name, :attribute_name) do |value, template, att|
                if att.raw.include?(prop_name.to_s)
                  raise ImportError.new("Value for #{prop_name} on attribute #{att.attribute_name} has duplicate definitions!")
                end

                att.send(:"#{prop_name}=", value)

                if att.attribute_type && att.link_model_name
                  if att.attribute_type == Etna::Clients::Magma::AttributeType::LINK && models.find_reciprocal(model_name: template.name, attribute: att).nil?
                    models.build_model(att.link_model_name).build_template.build_attributes.build_attribute(template.name).tap do |rec_att|
                      rec_att.attribute_name = rec_att.name = template.name
                      rec_att.display_name = prettify(template.name)
                      rec_att.description = prettify(template.name)
                      rec_att.attribute_type = Etna::Clients::Magma::AttributeType::COLLECTION
                      rec_att.link_model_name = template.name
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end
