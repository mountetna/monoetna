require 'json'
require 'openssl'
require 'net/http'
require 'singleton'

module Redcap
  class Client
    def initialize(host, token)
      @host = host
      @token = token
    end

    def get_data(request)
      request = request.merge(
        token: @token,
        format: 'json',
        returnFormat: 'json'
      )
      uri = URI("#{@host}/api/")
      response = Net::HTTP.post_form(uri, request)

      return nil unless response.content_type =~ %r!application/json!

      JSON.parse(response.body, symbolize_names: true)
    end

    def get_records(opts={})
      get_data(
        {
        content: 'record',
        rawOrLabel: 'label',
        rawOrLabelHeaders: 'raw',
        exportCheckboxLabel: 'true',
        exportSurveyFields: 'true',
        exportDataAccessGroups: 'true',
        }.merge(opts)
      )
    end

    def get_record_eavs(opts={})
      get_records({ type: 'eav' }.merge(opts))
    end

    def get_record_flat(opts={})
      get_records({ type: 'flat' }.merge(opts))
    end

    def records(form, events=false)
      eavs = get_record_eavs(form.fields)

      return nil unless eavs

      flat_records = get_record_flat(form.fields).map do |record|
        [ record[:record_id], record ]
      end.to_h

      records = eavs.group_by do |eav|
        events ?
          [ eav[:record], eav[:redcap_event_name] ]  :
          eav[:record]
      end.map do |record_id, record_eavs|
        [
          record_id,
          EavSet.new(
            record_eavs,
            form.form_name,
            flat_records[record_id],
            form.labels
          ).record
        ]
      end.to_h.compact

      records
    end

    class EavSet
      def initialize eavs, form_name, flat_record, labels
        @eavs = eavs
        @form_name = form_name
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
        eav[:field_name] != "#{@form_name}_complete" &&
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
