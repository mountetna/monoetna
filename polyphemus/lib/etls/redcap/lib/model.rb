module Redcap
  class Model
    def self.create(model_name, scripts, magma_template, redcap_template, salt)
      class_name = model_name.to_s.split('_').map(&:capitalize).join
      model_class = Kernel.const_defined?(class_name) ?  Kernel.const_get(class_name) : nil

      raise "No model class for #{model_name}" unless model_class

      model_class.new(scripts, magma_template, redcap_template, salt)
    end

    def initialize(scripts, magma_template, redcap_template, salt)
      @scripts = scripts.map do |script|
        Redcap::Script.new(self, script, redcap_template)
      end
      @magma_template = magma_template
      @salt = salt
      @offset_days = {}
    end

    def load(project)
      project.logger.write("Attempting to load model #{name}.\n")

      records = {}
      @scripts.each do |script|
        records.update(script.load(project))
      end

      records
    end

    def events?
      false
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

    def each_attribute(form)
      @terms.select{|t| t[:crf_form] == form }.each do |t|
        yield t[:attribute_name], t[:type], t[:crf_variable].to_sym
      end
    end

    def forms
      @terms.map{|l| l[:crf_form] }.compact.uniq.map
    end

    def name
      self.class.name.gsub(/(?<!^)[A-Z]/) do "_#$&" end.downcase.gsub("kernel::_", "")
    end

    def has_attribute?(att_name)
      @magma_template.attributes.attribute_keys.include?(att_name.to_s)
    end

    def attribute(att_name)
      @magma_template.attributes.attribute(att_name.to_s)
    end

    def cast_type(value, att_name, id)
      return nil unless value
      case attribute(att_name).attribute_type
      when "date_time"
        # eventually, we hope, magma will do this
        return nil if value.empty?

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
end
