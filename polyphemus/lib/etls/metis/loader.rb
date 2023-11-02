module Metis
  class Loader
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
