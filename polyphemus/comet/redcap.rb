require 'json'
require 'openssl'
require 'net/http'

class RedcapClient
  def initialize(token, host, salt)
    @token = token
    @host = host
    @salt = salt
    @offset_days = {}
  end

  def form_fields(form)
    ([ 'record_id' ] + template.select{|f| f[:form_name] == form}.map{|f| f[:field_name]}).map.with_index do |f,i|
      [ "fields[#{i}]", f ]
    end.to_h
  end

  def get_data(request)
    request = request.merge(
      token: @token,
      format: 'json',
      returnFormat: 'json'
    )
    uri = URI("https://#{@host}/api/")
    response = Net::HTTP.post_form(uri, request)

    return nil unless response.content_type =~ %r!application/json!

    JSON.parse(response.body, symbolize_names: true)
  end

  def template
    @template ||= get_data(content: 'metadata')
  end

  def forms
    @forms ||= template.map do |field|
      field[:form_name]
    end.uniq
  end

  def records(form, events=false)
    eavs = get_data(
      {
      content: 'record',
      type: 'eav',
      rawOrLabel: 'label',
      rawOrLabelHeaders: 'raw',
      exportCheckboxLabel: 'true',
      exportSurveyFields: 'true',
      exportDataAccessGroups: 'true',
      }.merge(form_fields(form))
    )

    return nil unless eavs

    records = eavs.group_by do |eav|
      events ?
        [ eav[:record], eav[:redcap_event_name] ]  :
        eav[:record]
    end.map do |record_id, record_eavs|
      [
        record_id,
        EavSet.new(
          record_eavs, form,
          labels
        ).record
      ]
    end.to_h.compact

    records
  end

  def clean_field(name)
    name.to_s.gsub(/_+[0-9]+$/,'').to_sym
  end

  def labels
    @labels ||= template.map do |t|
      [ t[:field_name].to_sym, t[:field_label] ]
    end.to_h
  end

  def offset_days(record_id)
    @offset_days[record_id] ||=
      begin
      # the offset in days is computed from hmac of the record_id
      signature = OpenSSL::HMAC.hexdigest(
        'SHA256',
        @salt,
        record_id
      )

      # we convert the hexadecimal string to a number in base 16.
      # A 64-character hex string becomes a 32 byte, 256 bit number
      # Divide by 2^256 to get a number between 0 and 1
      signature_fraction = signature.to_i(16).to_f/(1<<256)

      # offset days are computed as a number from 0 to 364
      (signature_fraction * 365).to_i
      end
  end

  class EavSet
    def initialize eavs, form, labels
      @eavs = eavs
      @form = form
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

    def make_record(record_eavs)
      record = record_eavs.reduce({}) do |rec,eav|
        if add_eav?(eav)
          rec[eav[:field_name].to_sym] = eav[:value]
        end
        rec
      end
      record.empty? ? nil : record
    end
  end
end
