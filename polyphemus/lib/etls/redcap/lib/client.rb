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
  end
end
