module Etna
  module Clients
    class Magma
      class ModelsCsv
        COLUMNS = [
            :comments,
            :model_name, :identifier, :parent_model_name, :parent_link_type,
            :attribute_name,
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
            options: [:validation, -> (s) { {"type" => "Array", "value" => s.split(',').map(&:chomp)} }],
            match: [:validation, -> (s) { {"type" => "Regexp", "value" => Regexp.new(s).source} }],
            attribute_group: :attribute_group,
            hidden: [:hidden, COLUMN_AS_BOOLEAN],
            unique: [:unique, COLUMN_AS_BOOLEAN],
        }

        def self.apply_csv_row(models = Models.new, row = {}, &err_block)

          models.tap do
            template = if (model_name = self.get_col_or_nil(row, :model_name))
              # Force the model to be most recently inserted key, since the state
              # of the current model is kept by key ordering.
              models.raw[model_name] = models.raw.delete(model_name)
              models.build_model(model_name).build_template.tap { |t| t.name = model_name }
            else
              last_model = models.raw.keys.last
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
            end
          end
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
            parent_model_attribute_by_model_name = parent_model.template.build_attributes.attribute(template.name)
            if parent_model_attribute_by_model_name && !reciprocal
              yield "Model #{template.parent} is linked as a parent to #{template.name}, but it already has an attribute named #{template.name} #{parent_model_attribute_by_model_name.raw}." if block_given?
            end

            parent_model.template.build_attributes.build_attribute(template.name).tap do |attr|
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

            if att.attribute_type == AttributeType::LINK && !models.find_reciprocal(model_name: template.name, attribute: att)
              models.build_model(att.link_model_name).build_template.build_attributes.build_attribute(template.name).tap do |rec_att|
                rec_att.attribute_name = template.name
                rec_att.display_name = self.prettify(template.name)
                rec_att.desc = self.prettify(template.name)
                rec_att.attribute_type = AttributeType::COLLECTION
                rec_att.link_model_name = template.name
              end
            end
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
          model_keys.each { |model_name| self.each_model_row(models, model_name, &block) }
        end

        def self.each_model_row(models, model_name, &block)
          return unless (model = models.model(model_name))
          # Empty link for better visual separation
          yield row_from_columns
          yield row_from_columns(model_name: model_name)

          unless model.template.parent.nil?
            return unless (parent_model = models.model(model.template.parent))
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

          model.template.attributes.all.each do |attribute|
            self.each_attribute_row(models, model, attribute, &block)
          end
        end

        def self.each_attribute_row(models = Models.new, model = Model.new, attribute = Attribute.new, &block)
          unless AttributeValidator.valid_add_row_attribute_types.include?(attribute.attribute_type) || attribute.attribute_type == AttributeType::IDENTIFIER
            return
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
              options: attribute.options&.join(', '),
              attribute_group: attribute.attribute_group,
              hidden: attribute.hidden,
              unique: attribute.unique,
          )
        end

        def self.row_from_columns(**columns)
          COLUMNS.map { |c| (columns[c] || '').to_s }
        end
      end
    end
  end
end
