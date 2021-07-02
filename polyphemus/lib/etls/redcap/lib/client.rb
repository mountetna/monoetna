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

      raise "Cannot fetch data from REDCap: #{response.body}." unless response.code == "200"

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

    class Record
      def initialize eavs, flat_record
        @eavs = eavs
        @flat_record = flat_record || {}
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
            rec[eav[:field_name].to_sym] = @flat_record[ eav[:field_name].to_sym ].yield_self do |flat_value|
              empty?(flat_value) ? eav[:value] : flat_value
            end
          end
          rec
        end
        record.empty? ? nil : record
      end
    end
  end
end
