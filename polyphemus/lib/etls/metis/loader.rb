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
        additionalProperties: {
          "$ref": "#/definitions/metis_model"
        }
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
              bucket_name: { type: "string" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              attribute_name: { type: "string" }
            },
            required: ["type", "bucket_name", "folder_path", "file_match", "attribute_name"]
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
              bucket_name: { type: "string" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              attribute_name: { type: "string" }
            },
            required: ["type", "bucket_name", "folder_path", "file_match", "attribute_name"]
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
              bucket_name: { type: "string" },
              folder_path: { type: "string" },
              file_match: { type: "string" },
              format: { enum: [ "tsv", "csv" ] },
              column_map: {
                type: "object",
                additionalProperties: { type: "string" }
              },
              extracted_columns: {
                type: "array",
                items: { type: "string" }
              }
            },
            required: ["type", "bucket_name", "folder_path", "file_match", "format" ]
          }
        }
      end
    end
  end
end
