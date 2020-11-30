require 'ostruct'

module Etna
  module Clients
    class Magma
      class ModelsCsv
        COPY_OPTIONS_SENTINEL = '$copy_from_original$'
        COLUMNS = [
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
        ]

        COLUMN_AS_BOOLEAN = -> (s) { ['true', 't', 'y', 'yes'].include?(s.downcase) }

        COLUMNS_TO_ATTRIBUTES = {
            attribute_name: :attribute_name,
            attribute_type: [:attribute_type, -> (s) { AttributeType.new(s) }],
            link_model_name: :link_model_name,
            description: :desc,
            display_name: :display_name,
            format_hint: :format_hint,
            restricted: [:restricted, COLUMN_AS_BOOLEAN],
            read_only: [:read_only, COLUMN_AS_BOOLEAN],
            options: [:validation, -> (s) { {"type" => "Array", "value" => s.split(',').map(&:strip)} }],
            match: [:validation, -> (s) { {"type" => "Regexp", "value" => Regexp.new(s).source} }],
            attribute_group: :attribute_group,
            hidden: [:hidden, COLUMN_AS_BOOLEAN],
            unique: [:unique, COLUMN_AS_BOOLEAN],
        }

        def self.apply_csv_row(changeset = ModelsChangeset.new, row = {}, &err_block)
          changeset.tap do
            models = changeset.models

            if (matrix_constant = self.get_col_or_nil(row, :matrix_constant))
              if matrix_constant.start_with?(COPY_OPTIONS_SENTINEL)
                matrix_constant = matrix_constant.slice((COPY_OPTIONS_SENTINEL.length)..-1)
                changeset.last_matrix_constant = matrix_constant
              else
                (changeset.matrix_constants[changeset.last_matrix_constant] ||= []) << matrix_constant
              end
            end

            template = if (model_name = self.get_col_or_nil(row, :model_name))
              changeset.last_model_key = model_name
              models.build_model(model_name).build_template.tap { |t| t.name = model_name }
            else
              last_model = changeset.last_model_key
              if last_model.nil?
                nil
              else
                models.model(last_model).build_template
              end
            end

            if (identifier = self.get_col_or_nil(row, :identifier))
              if template.nil?
                yield "Found identifier #{identifier} but no model_name had been given!"
                next
              end
              template.identifier = identifier
            end

            if (parent_model_name = self.get_col_or_nil(row, :parent_model_name))
              self.process_parent_model_name(template: template, parent_model_name: parent_model_name, models: models, &err_block)
            end

            if (parent_link_type = self.get_col_or_nil(row, :parent_link_type))
              self.process_parent_link_type(template: template, parent_link_type: parent_link_type, models: models, &err_block)
            end

            if (attribute_name = self.get_col_or_nil(row, :attribute_name))
              self.process_attribute(template: template, attribute_name: attribute_name, models: models, row: row, &err_block)
              if (new_attribute_name = self.get_col_or_nil(row, :new_attribute_name))
                self.process_new_attribute_name(template: template, changeset: changeset, new_attribute_name: new_attribute_name, attribute_name: attribute_name, &err_block)
              end
            end
          end
        end

        def self.process_new_attribute_name(template:, changeset:, new_attribute_name:, attribute_name:, &err_block)
          renames = changeset.build_renames(template.name)
          if renames.include?(attribute_name) && renames[attribute_name] != new_attribute_name
            if block_given?
              yield "Found multiple new_attribute_name values for #{template.name}'s #{attribute_name}': #{new_attribute_name} or #{renames[attribute_name]}?"
            end
          end

          renames[attribute_name] = new_attribute_name
        end

        def self.process_parent_model_name(template:, parent_model_name:, models:, &err_block)
          if template.nil?
            yield "Found parent_model_name #{parent_model_name} but no model_name had been given!" if block_given?
            return
          end

          if template.parent && !template.parent.empty? && template.parent != parent_model_name
            yield "Model #{template.name} was provided multiple parent_model_names: #{template.parent} and #{parent_model_name}" if block_given?
          end

          template.parent = parent_model_name

          template.build_attributes.build_attribute(template.parent).tap do |parent_att|
            parent_att.name = parent_att.attribute_name = parent_model_name
            parent_att.attribute_type = AttributeType::PARENT
            parent_att.link_model_name = parent_model_name
            parent_att.desc = self.prettify(parent_model_name)
            parent_att.display_name = self.prettify(parent_model_name)
          end
        end

        def self.process_parent_link_type(template:, parent_link_type:, models:, &err_block)
          if template.nil?
            yield "Found parent_link_type #{parent_link_type} but no model_name had been given!" if block_given?
            return
          end

          if template.parent.nil?
            yield "Found parent_link_type #{parent_link_type} but no parent_model_name had been given!" if block_given?
            return
          end

          reciprocal = models.find_reciprocal(model_name: template.name, link_attribute_name: template.parent)
          if reciprocal && reciprocal.attribute_type.to_s != parent_link_type
            yield "Model #{template.name} was provided multiple parent_link_types: #{reciprocal.attribute_type} and #{parent_link_type}" if block_given?
          end

          if reciprocal && reciprocal.attribute_name != template.name
            yield "Model #{template.name} is linked to #{template.parent}, but the reciprocal link is misnamed as '#{reciprocal.attribute_name}'." if block_given?
          end

          models.build_model(template.parent).tap do |parent_model|
            parent_model_attribute_by_model_name = parent_model.build_template.build_attributes.attribute(template.name)
            if parent_model_attribute_by_model_name && !reciprocal
              yield "Model #{template.parent} is linked as a parent to #{template.name}, but it already has an attribute named #{template.name} #{parent_model_attribute_by_model_name.raw}." if block_given?
            end

            parent_model.build_template.build_attributes.build_attribute(template.name).tap do |attr|
              attr.attribute_name = attr.name = template.name
              attr.attribute_type = parent_link_type
              attr.link_model_name = template.name
              attr.desc = self.prettify(template.name)
              attr.display_name = self.prettify(template.name)
            end
          end
        end

        def self.process_attribute(template:, attribute_name:, models:, row:, &err_block)
          if template.nil?
            yield "Found attribute #{attribute_name} but no model_name had been given!" if block_given?
            return
          end

          attributes = template.build_attributes
          existing_attribute = attributes.attribute(attribute_name)
          if existing_attribute
            if existing_attribute.attribute_type != AttributeType::COLLECTION
              yield "Attribute #{attribute_name} of model #{template.name} has duplicate definitions!" if block_given?
            end

            attributes.raw.delete(attribute_name)
          end

          attributes.build_attribute(attribute_name).tap do |att|
            COLUMNS_TO_ATTRIBUTES.each do |column, processor|
              next unless (value = self.get_col_or_nil(row, column))
              if processor.is_a?(Array)
                processor, f = processor
                value = f.call(value)
              end

              if !att.send(processor).nil? && !att.send(processor).empty?
                yield "Value for #{processor} on attribute #{attribute_name} has duplicate definitions!" if block_given?
              end

              att.send(:"#{processor}=", value)
            end

            if att.attribute_type == AttributeType::LINK && models.find_reciprocal(model_name: template.name, attribute: att).nil?
              models.build_model(att.link_model_name).build_template.build_attributes.build_attribute(template.name).tap do |rec_att|
                rec_att.attribute_name = template.name
                rec_att.display_name = self.prettify(template.name)
                rec_att.desc = self.prettify(template.name)
                rec_att.attribute_type = AttributeType::COLLECTION
                rec_att.link_model_name = template.name
              end
            end

            att.set_field_defaults!
          end
        end

        def self.prettify(name)
          name.split('_').map(&:capitalize).join(' ')
        end

        def self.get_col_or_nil(row, col)
          c = row[col]&.chomp
          return nil if c&.empty?
          c
        end

        def self.each_csv_row(models = Models.new, model_keys = models.model_keys.sort, &block)
          yield COLUMNS.map(&:to_s)

          self.ensure_parents(models, model_keys, &block)
          matrix_constants = {}
          model_keys.each { |model_name| self.each_model_row(models, model_name, matrix_constants, &block) }

          matrix_constants.each do |digest, options|
            yield row_from_columns
            yield row_from_columns(matrix_constant: COPY_OPTIONS_SENTINEL + digest)
            options.each do |option|
              yield row_from_columns(matrix_constant: option)
            end
          end
        end

        def self.ensure_parents(models, model_keys, &block)
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
              self.each_model_trunk_row(models, model_key, &block)
            end

            q.push(*models.model(model_key).template.all_linked_model_names)
          end
        end

        def self.each_model_trunk_row(models, model_name, &block)
          return unless (model = models.model(model_name))

          # Empty link for better visual separation
          yield row_from_columns
          yield row_from_columns(model_name: model_name)

          unless model.template.parent.nil?
            parent_model = models.model(model.template.parent)
            reciprocal = models.find_reciprocal(model: model, link_model: parent_model)

            yield row_from_columns(
                identifier: model.template.identifier,
                parent_model_name: model.template.parent,
                parent_link_type: reciprocal.attribute_type.to_s
            )
          else
            yield row_from_columns(
                identifier: model.template.identifier,
            )
          end
        end

        def self.each_model_row(models, model_name, matrix_constants, &block)
          return unless (model = models.model(model_name))

          self.each_model_trunk_row(models, model_name, &block)
          model.template.attributes.all.each do |attribute|
            self.each_attribute_row(models, model, attribute, matrix_constants, &block)
          end
        end

        def self.each_attribute_row(models = Models.new, model = Model.new, attribute = Attribute.new, matrix_constants = {}, &block)
          if attribute.attribute_type == AttributeType::IDENTIFIER
            # Identifiers for models whose parent link type ends up being a table are non configurable, so we don't
            # want to include them in the CSV.
            if models.find_reciprocal(model: model, link_attribute_name: model.template.parent)&.attribute_type == AttributeType::TABLE
              return
            end
          else
            return unless AttributeValidator.valid_add_row_attribute_types.include?(attribute.attribute_type)
          end

          options = attribute.options&.join(', ')
          if attribute.attribute_type == AttributeType::MATRIX
            # Matrix attribute validations are massive, and functional.  For now, we don't support showing and editing
            # them inline with this spreadsheet.  In the future, I think we should possibly introduce the concept of
            # CONSTANTS or Matrix Types that are managed separately.
            options = options || ''
            digest = Digest::MD5.hexdigest(options)
            matrix_constants[digest] ||= COLUMNS_TO_ATTRIBUTES[:options][1].call(options)["value"]

            options = COPY_OPTIONS_SENTINEL + digest
          end

          yield row_from_columns(
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

        def self.row_from_columns(**columns)
          COLUMNS.map { |c| (columns[c] || '').to_s }
        end

        class ModelsChangeset < Struct.new(:models, :renames, :matrix_constants, :last_matrix_constant, :last_model_key, keyword_init: true)
          def initialize(*args)
            super

            self.models ||= Models.new
            self.renames ||= {}
            self.matrix_constants ||= {}
          end
        end
      end
    end
  end
end
