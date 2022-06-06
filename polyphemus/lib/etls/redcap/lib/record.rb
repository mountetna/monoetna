module Redcap
  class Record
    def initialize eavs, flat_record, identifier_fields_data
      @eavs = eavs
      @flat_record = flat_record || {}
      @identifier_fields_data = identifier_fields_data || {}
    end

    def record
      make_record(@eavs)
    end

    private

    def empty?(value)
      value.nil? || value.empty?
    end

    def field_name(eav)
      eav[:field_name].to_s.gsub(/_?+[0-9]+$/,'').to_sym
    end

    def make_record(record_eavs)
      record = record_eavs.reduce({}) do |rec,eav|
        if !empty?(eav[:value])
          (rec[eav[:field_name].to_sym] ||= []).push(
            @flat_record[ eav[:field_name].to_sym ].yield_self do |flat_value|
              empty?(flat_value) ? eav[:value] : flat_value
            end
          )
        end
        rec
      end
      record.empty? ? nil : record.update(lift(@identifier_fields_data))
    end

    def lift(data)
      return {} if data.nil?

      # Make sure the values are wrapped in arrays, to match how loader behaves
      data.map do |key, value|
        [key, value.is_a?(Array) ? value : [value].compact]
      end.to_h
    end

  end
end
