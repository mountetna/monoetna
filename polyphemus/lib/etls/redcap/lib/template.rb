module Redcap
  class Template
    def initialize(template)
      @template = template
    end

    def forms
      @forms ||= @template.map do |field|
        field[:form_name]
      end.uniq.map(&:to_sym)
    end

    def field_names
      fields.keys
    end

    def fields
      @fields ||= @template.map do |f|
        [ f[:field_name], f ]
      end.to_h
    end

    def field_label(field_name)
      fields.dig(field_name,:field_label)
    end

    def field_note(field_name)
      fields.dig(field_name,:field_note)
    end
  end
end
