module Redcap
  class Payload
    attr_reader :eavs
    def initialize(eavs)
      @eavs = eavs

      @flat_records = flat_records.map do |record|
        [ record[:record_id], record ]
      end.to_h
    end

    def records
      @eavs.group_by do |eav|
        events ?
          [ eav[:record], eav[:redcap_event_name] ]  :
          eav[:record]
      end.map do |record_id, record_eavs|
        [
          record_id,
          EavSet.new(
            record_eavs, form,
            flat_records[record_id],
            labels
          ).record
        ]
      end.to_h.compact
    end

    class EavSet
      def initialize eavs, form, flat_record, labels
        @eavs = eavs
        @form = form
        @flat_record = flat_record || {}
        @labels = labels
      end

      def record
        if @eavs.first[:redcap_repeat_instance]
          records = @eavs.group_by do |eav|
            eav[:redcap_repeat_instance].to_s
          end.map do |repeat_id, record_eavs|
            next if empty?(repeat_id)
            make_record(record_eavs, false)
          end.compact
          records.empty? ? nil : records
        else
          make_record(@eavs)
        end
      end

      private

      def empty?(value)
        value.nil? || value == ''
      end

      def add_eav?(eav)
        eav[:field_name] != "#{@form}_complete" &&
          eav[:value] &&
          eav[:value] != '' &&
          ![ 'Not Available', 'Not Reported', 'Not reported' ].include?(eav[:value])
      end

      def field_name(eav)
        eav[:field_name].to_s.gsub(/_?+[0-9]+$/,'').to_sym
      end

      def make_record(record_eavs, use_flat=true)
        record = record_eavs.reduce({}) do |rec,eav|
          if add_eav?(eav)
            rec[eav[:field_name].to_sym] = best_value(
              eav, use_flat
            )
          end
          rec
        end
        record.empty? ? nil : record
      end

      def best_value(eav, use_flat)
        value = @flat_record[ eav[:field_name].to_sym ]

        !use_flat || empty?(value) ? eav[:value] : value
      end
    end
  end
end
