module Redcap
  class Script
    def initialize(model, script, template)
      @model = model
      @project = model.project
      @script = script
      @template = template
      @attributes ||= @script[:attributes].map do |att_name, att_value|
        [ att_name, Redcap::Value.new(att_value, template) ]
      end.to_h
    end

    def fields
      @fields ||= @attributes.values.map(&:field_name).uniq
    end

    def records
      eavs = @project.eav_records.select do | eav|
        fields.include?(eav[:field_name])
      end

      return {} unless eavs

      records = eavs.group_by do |eav|
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

      return records
    end

    def load
      update = {}

      records.each do |record_name, redcap_record|
        record_name, event_name = record_name if @model.events?
        next if record_name == "test"

        id = @model.identifier(record_name, event_name)

        next unless id

        next unless @attributes.values.all?{|v| v.valid?(redcap_record) }

        update[id] ||= {}

        @attributes.each do |att_name, att_value|
          next unless @model.has_attribute?(att_name)
          update[id][ att_name ] = @model.cast_type(
            att_value.to_value(redcap_record, id),
            att_name, id
          )
        end
      end

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
