module Metis
  class Loader
    class Error < Exception
    end
    class Script
      def initialize(model_config, script)
        @raw_script = script
        @model_config = model_config
      end

      def type
        @raw_script[:type]
      end

      def attribute_name
        @raw_script[:attribute_name]
      end

      def folder_path
        @raw_script[:folder_path]
      end

      def file_match
        @raw_script[:file_match]
      end

      def add_to_update(update, tail, metis=nil)
        case type
        when 'file'
          file_update(update, tail)
        when 'file_collection'
          file_collection_update(update, tail)
        when 'data_frame'
          data_frame_update(update, tail, metis)
        else
          raise Metis::Loader::Error.new("Invalid type for script #{type} for model #{@model_name}")
        end
      end

      def validate_rule
        if !@model_config.rule
          raise Metis::Loader::Error.new(
            "Cannot filter by file without a rule for #{@model_config.model_name}"
          )
        end

        @rule_match ||= Regexp.new(@model_config.rule[1..-2])
      end

      def file_update(update, tail)
        validate_rule

        files(tail).each do |file|
          name = identifier(file.file_path)
          next unless name
          update.update_revision(
            @model_config.model_name,
            name,
            {
              attribute_name => file.as_magma_file_attribute
            }
          )
        end
      end

      def identifier(path)
        match = @rule_match.match(path)

        return match ? match[0] : nil
      end

      def file_collection_update(update, tail)
        validate_rule

        grouped_files = files(tail).group_by do |file|
          identifier(file.file_path)
        end

        grouped_files.each do |name, name_files|
          next unless name
          update.update_revision(
            @model_config.model_name,
            name,
            {
              attribute_name => name_files.map(&:as_magma_file_attribute)
            }
          )
        end
      end

      def files(tail)
        tail.files.select do |file|
          ::File.fnmatch?(
            folder_path + '/**/' + file_match,
            file.file_path,
            ::File::FNM_PATHNAME | ::File::FNM_EXTGLOB
          )
        end
      end
    end

    class ModelConfig
      attr_reader :model_name
      attr_reader :model_def
      attr_reader :rule
      def initialize(model_name, model_config, model_def, rule)
        @model_name = model_name
        @raw_model_config = model_config
        @model_def = model_def
        @rule = rule
      end

      def scripts
        @scripts ||= (@raw_model_config[:scripts] || []).map do |script|
          Script.new(self, script)
        end
      end
    end

    class Config
      def initialize(config, params, model_defs, rules)
        @raw_config = config
        @params = params
        @model_defs = model_defs
        @rules = rules
      end

      def config
        @raw_config[:config] || {}
      end

      def models
        @models ||= (config[:models] || {}).map do |model_name, model_config|
          [ model_name, ModelConfig.new(
              model_name, model_config, @model_defs[model_name], @rules[model_name]
          ) ]
        end.to_h
      end

      def dry_run?
        ! @params[:commit]
      end

      def autolink?
        !!config[:autolink]
      end

      def project_name
        @raw_config[:project_name]
      end
    end

    def initialize(config = {}, rules = {}, params = {}, model_defs={})
      @config = Config.new(
        JSON.parse(JSON[config], symbolize_names: true),
        JSON.parse(JSON[params], symbolize_names: true),
        model_defs,
        rules
      )
    end

    def update_for(tail, metis_client=nil, models=nil)
      update = Etna::Clients::Magma::UpdateRequest.new(
        project_name: @config.project_name,
        dry_run: @config.dry_run?,
        autolink: @config.autolink?
      )

      @config.models.each do |model_name, model_config|
        model_config.scripts.each do |script|
          script.add_to_update(update, tail, metis_client)
        end
      end

      return update
    end

    def self.to_schema
      {
        "$schema": "http://json-schema.org/draft-07/schema#",
	title: "Metis Loader",
	description: "This loader takes data from Metis and imports it into Magma according to a specified mapping template.",
        definitions: [
          Metis::Loader::Model,
          Metis::Loader::File,
          Metis::Loader::FileCollection,
          Metis::Loader::DataFrame,
        ].map(&:to_schema).reduce(&:merge),
	type: "object",
        properties: {
          bucket_name: { type: "string" },
          autolink: { type: "boolean" },
          models: {
            type: "object",
            additionalProperties: {
              "$ref": "#/definitions/metis_model"
            }
          }
        },
        required: [ "bucket_name" ]
      }
    end

    class Model
      def self.to_schema
        {
          metis_model: {
            type: "object",
            properties: {
              scripts: {
                type: "array",
                items: {
                  oneOf: [
                    { "$ref": "#/definitions/metis_file" },
                    { "$ref": "#/definitions/metis_file_collection" },
                    { "$ref": "#/definitions/metis_data_frame" }
                  ]
                }
              }
            },
            additionalProperties: false
          }
        }
      end
    end
    class File
      def self.to_schema
        {
          metis_file: {
            type: "object",
            properties: {
              type: { const: "file" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              attribute_name: { type: "string" }
            },
            required: ["type", "folder_path", "file_match", "attribute_name"]
          }
        }
      end
    end
    class FileCollection
      def self.to_schema
        {
          metis_file_collection: {
            type: "object",
            properties: {
              type: { const: "file_collection" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              attribute_name: { type: "string" }
            },
            required: ["type", "folder_path", "file_match", "attribute_name"]
          }
        }
      end
    end
    class DataFrame
      def self.to_schema
        {
          metis_data_frame: {
            type: "object",
            properties: {
              type: { const: "data_frame" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              format: { enum: [ "tsv", "csv", "auto-detect" ] },
              blank_table: { type: "boolean" },
              column_map: {
                type: "object",
                minProperties: 2,
                additionalProperties: { type: "string" }
              },
              values_to_ignore: {type: "string"}
            },
            required: ["type", "folder_path", "file_match", "format", "column_map"]
          }
        }
      end
    end
  end
end
