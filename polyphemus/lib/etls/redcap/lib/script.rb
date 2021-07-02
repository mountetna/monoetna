module Redcap
  class Script
    def initialize(model, script, template)
      @model = model
      @project = model.project
      @template = template
      @attributes ||= script[:attributes].map do |att_name, att_value|
        [ att_name, Redcap::Value.new(att_name, att_value, template) ]
      end.to_h
      @each_entities = script[:each]
    end

    def fields
      @fields ||= @attributes.values.map(&:field_name).uniq
    end

    def group_by_iter(eavs)
      eavs.group_by do |eav|
        (@each_entities || @model.each_entities).map do |ent|
          case ent
          when :record
            eav[:record]
          when :event
            eav[:redcap_event_name]
          when :repeat
            [ eav[:redcap_repeat_instrument], eav[:redcap_repeat_instance] ]
          when :value
            eav[:value]
          end
        end.compact
      end
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

      @model.existing_records.each do |magma_record_name, magma_record|
        next unless @attributes.values.all?{|v| v.valid?(magma_record) }

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
        magma_record_name = @model.identifier(*record_id)

        next unless magma_record_name

        next unless @attributes.values.all?{|v| v.valid?(redcap_record) }

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
  end
end
