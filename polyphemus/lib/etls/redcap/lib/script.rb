module Redcap
  class Script
    def initialize(model, script)
      @model = model
      @script = script
    end

    def forms
      @forms ||= @script[:forms].map do |form_name, form|
        [ form_name, Redcap::Form.new(@model, form_name, form) ]
      end.to_h
    end

    def load(client)
      records = {}

      forms.each do |form_name, form|
        raise "Missing form #{form_name}" unless client.has_form?(form_name)

        puts "Processing form #{form_name}"

        form.load(client, records)
      end

      puts "Patching unfilled attributes"

      records.each do |id, record|
        record.compact!

        @model.patch(id, record) unless record.empty?
      end
        
      records.select do |id,record|
        !record.empty?
      end.to_h
    end
  end
end
