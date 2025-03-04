module Metis
  class Loader
    class Error < Exception
    end
    class Script
      def initialize(config, model_name, script)
        @config = config
        @raw_script = script
        @model_name = model_name
      end

      def self.script_params(*params)
        params.each do |param|
          self.define_method param do
            @raw_script[param]
          end
        end
      end

      script_params :type, :attribute_name, :folder_path,
        :file_match, :column_map, :format, :values_to_ignore,
        :blank_table

      def model
        @config.model_defs.model(model_name)
      end

      def model_name
        @model_name
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
          raise Metis::Loader::Error.new("Invalid type for script #{type} for model #{model_name}")
        end
      end

      def is_table?
        @config.model_defs.is_table?(model_name)
      end

      def file_update(update, tail)
        files(tail).each do |file|
          name = @config.identifier(model_name, file.file_path)
          next unless name
          update.update_revision(
            model_name,
            name,
            {
              attribute_name => file.as_magma_file_attribute
            }
          )
        end
      end

      def table_separator(file)
        case format
        when 'csv'
          return ','
        when 'tsv'
          return "\t"
        when 'auto-detect'
          if file.file_path =~ /csv$/
             return ','
          elsif file.file_path =~ /tsv$/
             return "\t"
          end
        end
      end

      def data_frame_update(update, tail, metis)
        model_attributes = model.template.attributes
        id = model.template.identifier

        if !is_table? && !column_map.has_key?(model.template.identifier.to_sym)
          raise Metis::Loader::Error.new(
            "Identifier attribute is missing from the 'column_map' of #{model_name} data_frame loader."
          )
        end

        if is_table? && !column_map.has_key?(model.template.parent.to_sym)
          raise Metis::Loader::Error.new(
            "Parent attribute is missing from the 'column_map' of #{model_name} data_frame loader."
          )
        end

        missing = (column_map.keys - model.template.attributes.attribute_keys.map(&:to_sym))

        if !missing.empty?
          raise Metis::Loader::Error.new(
            "'column_map' of #{model_name} data_frame loader targets attribute(s) that don't exist: #{missing.join(', ')}."
          )
        end

        # Parse Files
        files(tail).each do |file|
          file_contents = String.new
          metis.download_file(file) do |chunk|
            file_contents << chunk
          end

          head, *rows = CSV.parse(StringIO.new(file_contents, 'r:bom|utf-8'), col_sep: table_separator(file))
          data = Daru::DataFrame.rows(rows, order: head.map(&:to_sym))

          if !data.any? || data.ncols < 2
            raise Metis::Loader::Error.new(
              "#{file.file_name} seems to have fewer than 2 columns. Check the 'format' configuration for this data_frame loader."
            )
          end

          missing = column_map.values - data.vectors.map(&:to_s)
          if !missing.empty?
            raise Metis::Loader::Error.new("#{file.file_name} is missing column(s) targetted by #{model_name} data_frame loader 'column_map': #{missing.join(', ')}.")
          end

          # not sure this is required any more
          #if data.any?(:row) { |r| r.any?(&:nil?) }
          #  raise Metis::Loader::Error.new("#{file.file_name} has unexpected empty values after all parsing. Data rows may be shorter than the column row indicates.")
          #end
          
          # Trim to mapped columns and convert to attribute names
          data.vectors = Daru::Index.new(
            data.vectors.map do |column_name|
              column_map.invert[ column_name.to_s ] || column_name
            end
          )

          data.vectors.to_a.each do |name|
            data.delete_vector name unless column_map.has_key?(name)
          end

          # Blank data equaling values_to_ignore by setting as nil
          if values_to_ignore
            data.replace_values(values_to_ignore.split(","), nil)
          end

          # Determine Updates
          if is_table?
            data[:__temp__] = data.index.map { |i| "::temp-id-#{i}" }
            data = data.set_index(:__temp__)

            if blank_table
              parent_model_name = model.template.parent.to_sym
              data[ parent_model_name ].each do |parent_name|
                next unless @config.identifier(parent_model_name, parent_name)
                update.update_revision(
                  parent_model_name.to_s,
                  parent_name,
                  {
                    model_name => data.where( data[parent_model_name].eq(parent_name) ).index.to_a
                  }
                )
              end
            end
            data.each_row_with_index do |attributes,name|
              update.update_revision(model_name.to_s, name, attributes.to_h.compact) 
            end
          else
            data = data.set_index(model.template.identifier.to_sym)
            data.each_row_with_index do |attributes,name|
              next if @config.rules[model_name] && !@config.identifier(model_name, name)
              update.update_revision(model_name.to_s, name, attributes.to_h.compact) 
            end
          end
        end
      end

      def file_collection_update(update, tail)
        grouped_files = files(tail).group_by do |file|
          @config.identifier(model_name, file.file_path)
        end

        grouped_files.each do |name, name_files|
          next unless name
          update.update_revision(
            model_name,
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

    class Config
      attr_reader :rules
      attr_reader :model_defs

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
          [
            model_name,
            (model_config[:scripts] || []).map do |script|
              Script.new(self, model_name, script)
            end
          ]
        end.to_h
      end

      def rule_match(model_name)
        @rule_match ||= {}
        @rule_match[model_name] ||= Regexp.new(@rules[model_name][1..-2])
      end

      def identifier(model_name, path)
        raise Metis::Loader::Error.new(
            "Cannot filter without a rule for #{model_name}"
        ) unless @rules[model_name]

        match = rule_match(model_name).match(path)

        return match ? match[0] : nil
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

    attr_reader :config

    def initialize(config = {}, rules = {}, params = {}, model_defs={})
      @config = Config.new(
        JSON.parse(JSON[config], symbolize_names: true),
        JSON.parse(JSON[params], symbolize_names: true),
        model_defs,
        rules
      )
    end

    def update_for(tail, metis_client=nil)
      update = Etna::Clients::Magma::UpdateRequest.new(
        project_name: @config.project_name,
        dry_run: @config.dry_run?,
        autolink: @config.autolink?
      )

      @config.models.each do |model_name, scripts|
        scripts.each do |script|
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
