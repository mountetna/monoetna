module Redcap
  class Form
    attr_reader :form_name

    def initialize(model, form_name, form)
      @model = model
      @form_name = form_name
      @form ||= form.map do |att_name, att_value|
        [ att_name, Redcap::Value.new(att_value) ]
      end.to_h
    end

    def load(client, records)
      client.records(@form_name, @model.events?).each do |record_name, redcap_record|
        record_name, event_name = record_name if @model.events?
        next if record_name == "test"

        id = @model.identifier(record_name, event_name)

        next unless id

        next unless @form.values.all?{|v| v.valid?(redcap_record) }

        records[id] ||= {}

        @form.each do |att_name, form_value|
          next unless @model.has_attribute?(att_name)
          records[id][ att_name ] = @model.cast_type(
            form_value.to_value(redcap_record, id, client.template),
            att_name, id
          )
        end
      end
    end
  end
end
