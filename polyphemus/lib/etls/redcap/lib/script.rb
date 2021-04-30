module Redcap
  class Script
    def initialize(model, script, template)
      @model = model
      @project = model.project
      @template = template
      @attributes ||= script.map do |att_name, att_value|
        [ att_name, Redcap::Value.new(att_name, att_value, template) ]
      end.to_h
    end

    def fields
      @fields ||= @attributes.values.map(&:field_name).uniq
    end

    def redcap_records
      return @redcap_records if @redcap_records

      eavs = @project.eav_records.select do | eav|
        fields.include?(eav[:field_name])
      end

      return {} unless eavs

      @redcap_records = eavs.group_by do |eav|
        @model.events? ?  [ eav[:record], eav[:redcap_event_name] ] : eav[:record]
      end.map do |record_id, record_eavs|
        [
          record_id,
          Redcap::Client::Record.new(
            record_eavs,
            @project.flat_records[record_id]
          ).record
        ]
      end.to_h.compact

      return @redcap_records
    end

    def inverse_load
      update = {}

      @model.existing_records.each do |record_name, magma_record|
        next unless @attributes.values.all?{|v| v.valid?(magma_record) }

        id = record_name

        redcap_record = redcap_records[ @model.redcap_id(*record_name) ]

        next unless redcap_record

        update[id] ||= {}

        update_record(update[id], id, redcap_record)
      end

      return patched(update)
    end

    def load
      update = {}

      redcap_records.each do |record_name, redcap_record|
        record_name, event_name = record_name if @model.events?
        next if record_name == "test"

        id = @model.identifier(record_name, event_name)

        next unless id

        next unless @attributes.values.all?{|v| v.valid?(redcap_record) }

        update[id] ||= {}

        update_record(update[id], id, redcap_record)
      end

      return patched(update)
    end

    def update_record(record, id, redcap_record)
      @attributes.each do |att_name, att_value|
        next unless @model.has_attribute?(att_name)

        next if att_value.none?

        record[ att_name ] = @model.cast_type(
          att_value.to_value(redcap_record, id),
          att_name, id
        )
      end
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
  end
end
