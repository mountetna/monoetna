module Redcap
  class Model
    def self.to_schema
      {
        redcap_model: {
          type: "object",
          properties: {
            each: { "$ref": "#/definitions/each" },
            invert: { type: "boolean" },
            identifier_fields: {
              type: "array",
              items: {
                type: "string"
              }
            },
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
      Kernel.const_set(model_name, Class.new(Redcap::Model))
    end

    def self.create(project, model_name, model_config, magma_template, redcap_template)
      class_name = model_name.to_s.split('_').map(&:capitalize).join
      model_class = Kernel.const_defined?(class_name) ?  Kernel.const_get(class_name) : Redcap::Model

      raise "No model class for #{model_name}" unless model_class

      model_class.new(model_name, project, model_config, magma_template, redcap_template)
    end

    attr_reader :scripts, :project, :model_name

    def initialize(model_name, project, config, magma_template, redcap_template)
      @model_name = model_name
      @project = project
      @config = config
      @scripts = config[:scripts].map do |script|
        Redcap::Script.new(self, script, redcap_template)
      end
      @magma_template = magma_template
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
          model_name: model_name,
          attribute_names: attribute_names,
          record_names: "all",
        )
      ).models.model(model_name.to_s).documents.raw
    end

    def attribute_names
      @config[:attributes] || 'identifier'
    end

    def identifier_fields
      @config[:identifier_fields] || []
    end

    def identifier(record_id, identifier_fields: nil)
      unless @magma_template.identifier.present?
        return [
            "::temp", *record_id, rand(36**8).to_s(36)
        ].compact.join('-')
      end

      identifier_fields.present? ? identifier_fields.values.first : record_id
    end

    def patch(id, record)
    end

    def ensure_containing_records?
      false
    end

    def ensure_records(current_records)
      {}
    end
  
    def load
      @project.logger.write("Attempting to load model #{model_name}.\n")

      records = {}
      @scripts.each do |script|
        script_records = invert? ? script.inverse_load : script.load

        script_records.each do |record_name, record|
          (records[record_name] ||= {}).update( record )
        end
      end

      records
    end

    def redcap_id(magma_record_name, magma_record)
      nil
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
      when "float"
        return value.to_f
      when "integer"
        return value.to_i
      when "boolean"
        yes_values = ["Yes", "True"]
        no_values = ["No", "False"]
        return yes_values.include?(value) ? true : no_values.include?(value) ? false : nil
      else
        return value
      end
    end
  end
end
