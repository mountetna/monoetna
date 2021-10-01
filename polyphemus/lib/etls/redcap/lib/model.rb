module Redcap
  class Model
    def self.to_schema
      {
        redcap_model: {
          type: "object",
          properties: {
            each: { "$ref": "#/definitions/each" },
            invert: { type: "boolean" },
            scripts: {
              type: "array",
              items: { "$ref": "#/definitions/script" }
            }
          },
          additionalProperties: false,
          required: [ "scripts" ]
        }
      }
    end

    def self.define(model_name, &block)
      return Kernel.const_get(model_name) if Kernel.const_defined?(model_name)

      # Set some default methods for each model
      Kernel.const_set(model_name, Class.new(Redcap::Model) {
        def identifier(*record_id)
          [
              "::temp", *record_id, rand(36**8).to_s(36)
          ].compact.join('-')
        end

        def patch(id, record)
        end
      })
    end

    def self.create(project, model_name, model_config, magma_template, redcap_template, salt)
      class_name = model_name.to_s.split('_').map(&:capitalize).join
      model_class = Kernel.const_defined?(class_name) ?  Kernel.const_get(class_name) : nil

      raise "No model class for #{model_name}" unless model_class

      model_class.new(model_name, project, model_config, magma_template, redcap_template, salt)
    end

    attr_reader :scripts, :project

    def initialize(model_name, project, config, magma_template, redcap_template, salt)
      @model_name = model_name
      @project = project
      @config = config
      @scripts = config[:scripts].map do |script|
        Redcap::Script.new(self, script, redcap_template)
      end
      @magma_template = magma_template
      @salt = salt
      @offset_days = {}
    end

    def invert?
      @config[:invert]
    end

    def each_entities
      @entities = @config[:each].map do |ent|
        Redcap::Entity.create(ent)
      end
    end

    def existing_records
      @existing_records ||= project.magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(
          project_name: @project.project_name,
          model_name: @model_name,
          attribute_names: attribute_names,
          record_names: "all",
        )
      ).models.model(@model_name.to_s).documents.raw
    end

    def attribute_names
      @config[:attributes] || 'identifier'
    end

    def load
      @project.logger.write("Attempting to load model #{name}.\n")

      records = {}
      @scripts.each do |script|
        records.update(
          invert? ? script.inverse_load : script.load
        )
      end

      records
    end

    def offset_id(record_id)
      record_id
    end

    def redcap_id(magma_record_name, magma_record)
      nil
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

        begin
          return (DateTime.parse(value) - offset_days(id)).iso8601[0..9]
        rescue ArgumentError
          return nil
        end
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
