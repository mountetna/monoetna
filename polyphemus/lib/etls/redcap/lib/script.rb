module Redcap
  class Script
    def initialize(model, script, redcap_template)
      @model = model
      @script = script
      @redcap_template = redcap_template
    end

    def forms
      @forms ||= @script[:forms].map do |form_name, form|
        [ form_name, Redcap::Form.new(@model, form_name, form, @redcap_template) ]
      end.to_h
    end

    def load(project)
      records = {}

      forms.each do |form_name, form|
        raise "Missing form #{form_name}" unless project.has_form?(form_name)

        project.logger.write("Processing form #{form_name}.\n")

        form.load(project, records)
      end

      project.logger.write("Patching unfilled attributes.\n")

      records.each do |id, record|
        record.compact! unless project.strict

        @model.patch(id, record) unless record.empty?
      end

      records.select do |id,record|
        !record.empty?
      end.to_h
    end
  end
end
