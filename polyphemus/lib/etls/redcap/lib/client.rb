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

      # For some reason :redcap_repeat_instance returns as int
      #   for non-eav, but str for eav when the value > 1...
      #   so we'll cast the values in this field to int so can do
      #   accurate mapping between eav and flat records.
      JSON.parse(response.body, symbolize_names: true).map do |record|
        record[:redcap_repeat_instance] = record[:redcap_repeat_instance].to_i unless (record[:redcap_repeat_instance].nil? || record[:redcap_repeat_instance] == '')

        record
      end
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
  end
end
