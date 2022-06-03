module Redcap
  class Script
    def self.to_schema
      {
        script: {
          type: "object",
          properties: {
            attributes: {
              type: "object",
              additionalProperties: {
                oneOf: [
                  { type: "string" },
                  { "$ref": "#/definitions/attribute_value" }
                ]
              }
            },
            each: { "$ref": "#/definitions/each" },
            filters: {
              type: "array",
              items: {
                "$ref": "#/definitions/filter_value"
              },
              uniqueItems: true
            }
          },
          additionalProperties: false,
          required: [ "attributes" ]
        }
      }
    end

    def initialize(model, script, template)
      @model = model
      @project = model.project
      @template = template
      @attributes ||= script[:attributes].map do |att_name, att_value|
        [ att_name, Redcap::Value.new(att_name, att_value, template) ]
      end.to_h
      @each_entities = script[:each]&.map do |ent|
        Redcap::Entity.create(ent)
      end
      @filter_entities = script[:filters]&.map do |filter|
        Redcap::Filter.create(filter)
      end
    end

    def fields
      @fields ||= @attributes.values.map(&:field_name).concat(@filter_entities&.map(&:field_name) || []).uniq.compact
    end

    def each_entities
      @each_entities || @model.each_entities
    end

    def group_by_iter(eavs)
      fields, groups = eavs.group_by do |eav|
        each_entities.map do |ent|
          ent.key(eav)
        end
      end.reject do |record_id, record_eavs|
        record_id.any?(&:nil?)
      end.partition do |record_id, record_eavs|
        record_id.any?{|i| i == :field }
      end.map(&:to_h)

      unless fields.empty?
        groups.each do |record_id, record_eavs|
          masked_record_id = record_id.map.with_index do |id, i|
            each_entities[i].is_a?(Redcap::Entity::Field) ? :field : id
          end
          record_eavs.concat(fields[masked_record_id]) if fields[masked_record_id]
        end
      end

      groups
    end

    def key_entities
      @key_entities ||= each_entities.select(&:flat_key?)
    end

    def flat_record(record_id)
      @flat_records ||= @project.flat_records.group_by do |record|
        key_entities.map{|e| e.key(record)}
      end

      # there may be several matching flat_records
      flat_record_for_id = @flat_records[
        key_entities.map.with_index {|e,i| e.flat_key? ? [ record_id[i] ] : nil }.compact.flatten(1)
      ]&.first

      # For repeating instruments, need to merge this with the
      #   base, non-repeating record to catch any specific
      #   identifiers...
      flat_record_for_id&.update(alternate_identifier_fields(record_id))
    end

    def redcap_records
      return @redcap_records if @redcap_records

      eavs = @project.eav_records.select do | eav|
        fields.include?(eav[:field_name])
      end

      return {} unless eavs

      @redcap_records = group_by_iter(eavs).map do |record_id, record_eavs|
        [
          record_id,
          Redcap::Record.new(
            record_eavs,
            flat_record(record_id),
            some(identifier_fields_data(record_id))
          ).record
        ]
      end.to_h.compact

      return @redcap_records
    end

    def inverse_load
      update = {}

      @model.existing_records.each do |magma_record_name, magma_record|

        next unless @attributes.values.all?{|v| v.valid_magma?(magma_record) }
        redcap_record = redcap_records[ @model.redcap_id(magma_record_name, magma_record) ]

        next unless redcap_record

        update.merge!(
          magma_record_name => update_record(magma_record_name, redcap_record)
        )
      end

      return patched(update)
    end

    def load
      update = {}

      redcap_records.each do |record_id, redcap_record|
        magma_record_name = @model.identifier(*record_id, identifier_fields: identifier_fields_data(record_id))

        next unless magma_record_name

        next unless @attributes.values.all?{|v| v.valid_redcap?(redcap_record) }

        if @filter_entities
          next unless @filter_entities.all? { |f| f.allow_redcap?(redcap_record) }
        end

        update.merge!(
          magma_record_name => update_record(magma_record_name, redcap_record)
        )
      end

      return patched(update)
    end

    def update_record(magma_record_name, redcap_record)
      @attributes.map do |att_name, att_value|
        next unless @model.has_attribute?(att_name)

        next if att_value.none?

        [
          att_name,
          @model.cast_type(
            att_value.to_value(redcap_record, magma_record_name),
            att_name, magma_record_name
          )
        ]
      end.compact.to_h
    end

    def patched(update)
      update.each do |id, record|
        record.compact! unless @project.strict

        @model.patch(id, record) unless record.empty?
      end

      update.select do |id,record|
        !record.empty?
      end.to_h
    end

    def identifier_fields_data(record_id)
      flat_record(record_id)&.slice(*identifier_fields)
    end

    def identifier_fields
      @identifier_fields ||= @model.identifier_fields.map do |field|
        field.to_sym
      end
    end

    def has_identifier_fields?
      @has_identifier_fields ||= @model.identifier_fields.length > 0
    end

    def is_repeating?
      @is_repeating ||= each_entities.any? { |e| e.is_a?(Redcap::Entity::Repeat) }
    end

    def alternate_identifier_fields(record_id)
      return {} unless is_repeating? && has_identifier_fields?

      return {} unless @flat_records
      
      non_repeating_key = non_repeating_record_key(record_id)
      
      return {} if non_repeating_key.nil?

      @flat_records[non_repeating_key]&.first.slice(*identifier_fields)
    end

    def some(data)
      return {} if data.nil?

      # Make sure the values are wrapped in arrays, to match how loader behaves
      data.map do |key, value|
        [key, value.is_a?(Array) ? value : [value].compact]
      end.to_h
    end

    def non_repeating_record_keys
      # flat record keys are Tuples in the form of
      # [ redcap_id, [ repeating_instrument, instance_id ] ]
      @non_repeating_record_keys ||= @flat_records.keys.select do |key|
        key.last.first == ""
      end 
    end

    def non_repeating_record_key(record_id)
      non_repeating_record_keys.select do |key|
        key.first == record_id.first
      end.first
    end
  end
end
