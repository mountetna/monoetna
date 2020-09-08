#!/usr/bin/env ruby

require './redcap'
require 'yaml'
require 'csv'
require 'etna'

class CometModel
  def self.create(model_name, terms, salt)
    model_class = Kernel.const_defined?(model_name.capitalize) ?  Kernel.const_get(model_name.capitalize) : nil

    raise "No model class for #{model_name}" unless model_class

    model_class.new(terms, salt)
  end

  def guess(att_name, att_type, crf_variable, id, redcap_record)
    nil
  end

  def offset_id(record_id)
    record_id
  end

  def offset_days(record_id)
    @offset_days[record_id] ||=
      begin
      # the offset in days is computed from hmac of the record_id
      signature = OpenSSL::HMAC.hexdigest(
        'SHA256',
        @salt,
        offset_id(record_id)
      )

      # we convert the hexadecimal string to a number in base 16.
      # A 64-character hex string becomes a 32 byte, 256 bit number
      # Divide by 2^256 to get a number between 0 and 1
      signature_fraction = signature.to_i(16).to_f/(1<<256)

      # offset days are computed as a number from 0 to 364
      (signature_fraction * 365).to_i
      end
  end

  def initialize(terms, salt)
    @terms = terms
    @salt = salt
    @offset_days = {}
  end

  def each_attribute(form)
    @terms.select{|t| t[:crf_form] == form }.each do |t|
      yield t[:attribute_name], t[:type], t[:crf_variable].to_sym
    end
  end

  def forms
    @terms.map{|l| l[:crf_form] }.compact.uniq
  end

  def events?
    false
  end

  def name
    self.class.name.gsub(/(?<!^)[A-Z]/) do "_#$&" end.downcase
  end

  def cast_type(value, type, id)
    return nil unless value
    case type
    when "date_time"
      # eventually, we hope, magma will do this
      return (DateTime.parse(value) - offset_days(id)).iso8601[0..9]
    when "float"
      return value.to_f
    when "integer"
      return value.to_i
    when "boolean"
      return value == "Yes" ? true : value == "No" ? false : nil
    else
      return value
    end
  end
end

class Patient < CometModel
  def identifier(record_name, event_name=nil)
    "MVIR1-HS#{record_name.to_i - 1000}"
  end

  CONSENT={
    "100" => 'Full Study',
    "90" => 'Current Samples/Data Only',
    "75" => 'Permanent Waiver',
    "50" => 'Initial Waiver',
    "0" => 'No Samples/Data'
  }
  def patch(id, record)
    return if record.empty?
    record["project"] = "COMET"
    record["consent"] = CONSENT[record["consent"]]
    record["age"] = [ 89, record["age"] ].min if record["age"]
  end
end


class Timepoint < CometModel
  def identifier(record_name, event_name)
    day = event_name[ /^D[0-9]+/ ]

    day ? "MVIR1-HS#{record_name.to_i - 1000}-#{day}" : nil
  end

  def events?
    true
  end

  def offset_id(record_id)
    record_id.sub(/-D[0-9]+$/,'')
  end

  def patch(id, record)
    return if record.empty?
    record["patient"] = id.split(/-/)[0..1].join('-')
    record["day"] = id[/-D([0-9]+)$/,1].to_i
    record["study_day"] = (record["day"] == 4) || (record["day"] % 7 == 0)
  end
end

class Treatment < CometModel
  def identifier(record_name, event_name)
    "::temp#{record_name}"
  end

  def each_variant(crf_variable)
    m = crf_variable.match(/\[(?<code>.*)\]/)
    if !m
      yield(crf_variable)
      return
    end

    m
  end

  def guess(att_name, att_type, crf_variable, id, redcap_record)
    return nil unless att_name == "name"

    each_variant(crf_variable) do crf_var
      require 'pry'
      binding.pry
    end
    nil
  end


  def patch(id, record)
    return if record.empty?
  end
end

config = YAML.load(File.read "config.yml")

r = RedcapClient.new(config[:token], config[:redcap_host], config[:dateshift_salt])

comet = CSV.read("comet.tsv", col_sep:"\t").yield_self do |t|
  header = t.shift.compact.map(&:to_sym)
  t.map do |r|
    row = header.zip(r.slice(0,header.length)).to_h
    row[:crf_form] = row[:crf_form]&.gsub(/\s/, '_')&.downcase
    row
  end
end

puts "Forms: #{r.forms}"

models = comet.map{|l| l[:model_name]}.compact.uniq

revisions = {}

models.each do |model_name|
  next unless model_name == "timepoint"

  model = CometModel.create(model_name, comet.select {|a| a[:model_name] == model_name }, config[:dateshift_salt])

  puts "Attempting to load model #{model.name}"

  missing_forms = model.forms - r.forms

  puts "Found #{model.forms.size-missing_forms.size} out of #{model.forms.size} required forms"

  puts "Missing #{missing_forms.join(', ')}" unless missing_forms.empty?

  puts "Creating records using each form"

  records = {}
  (model.forms & r.forms).each do |form|
    puts "Processing form #{form}"
    r.records(form, model.events?).each do |record_name, redcap_record|
      record_name, event_name = record_name if model.events?
      next if record_name == "test"

      id = model.identifier(record_name, event_name)
      next unless id
      records[id] ||= {}

      model.each_attribute(form) do |att_name, att_type, crf_variable|
        records[id][ att_name ] = model.guess(att_name, att_type, crf_variable, id, redcap_record) || model.cast_type(
          redcap_record[ crf_variable ], att_type, id
        )
      end
    end
  end

  puts "Patching unfilled attributes"

  records.each do |id, record|
    record.compact!

    model.patch(id, record)
  end
    
  revisions[model.name] = records.select do |id,record|
    !record.empty?
  end.to_h
end

puts "Posting revisions"

e = Etna::Client.new('https://magma.ucsf.edu', ENV['TOKEN'])

response = e.update(project_name: "mvir1", revisions: revisions)
