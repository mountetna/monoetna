require 'yaml'
require_relative 'contexts'
require_relative 'ruby'

module Etna
  module Codegen
    class Swagger < Project
      def initialize(logger = Etna::Logger.new('/proc/self/fd/2', 0, 1024))
        super('etna')
        @logger = logger
        @swagger_cache = {}
      end

      def load_swagger(path)
        @swagger_cache ||= SwaggerFile.new(self, path)
      end

      def generate_api(swagger_file)
        etna_client(app_name)
      end

      def generate_all
        Dir[File.expand_path("../../swagger/*.yml", __FILE__)].each do |api_file|
          generate_api(load_swagger(api_file))
        end
      end

      private

      def etna_client(swagger_file)
        app_name = swagger_file.basename(".*.yml")
        build_file("clients/#{app_name}/client.rb") do |key|
          Ruby::EtnaClientFile.new(key.update(app_name: app_name, swagger_file: swagger_file, scopes: {module: "Etna::Clients"}))
        end
      end
    end

    class SwaggerFile
      def initialize(swagger, abs_path)
        @swagger = swagger
        @abs_path = abs_path
      end

      def basename(ext)
        File.basename(@abs_path, ext)
      end

      def contents
        @contents ||= begin
          YAML.load_file(@abs_path)
        end
      end

      def load_file(path)
        @swagger.load_file(path)
      end
    end
  end
end
