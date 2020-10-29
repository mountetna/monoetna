# Creates a project from a JSON file.
# This workflow:
# 1) Creates the project in Janus.
# 2) Adds the user as a project administrator to Janus.
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in Magma.
# 5) Iterates through the models and attributes and creates those in Magma.

require 'base64'
require 'json'
require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative '../../janus/models'
require_relative './json_validators'
require_relative './json_converters'

module Etna
  module Clients
    class Magma
      class CreateProjectModelsWorkflow < Struct.new(:magma_client, keyword_init: true)
        attr_reader :project, :converter

        # def initialize(**params)
        #   user_json = JSON.parse(File.read(filepath))
        #
        #   magma_json = Etna::Clients::Magma::ProjectConverter.convert_project_user_json_to_magma_json(user_json)
        #   @project = Etna::Clients::Magma::Project.new(magma_json)
        #
        #   @validator = Etna::Clients::Magma::ProjectValidator.new(project)
        #   @validator.validate
        #
        #   raise "Project JSON has errors: #{@validator.errors}" unless @validator.valid?
        #
        #   @converter = ProjectConverter.new(project)
        #   @converter.convert!
        # end

        def process_file(starting_model_tsv) end

        def each_tsv_row_for_template(project_name, target_model = 'project', &block)
          models = magma_client.retrieve(RetrievalRequest.new(project_name: project_name, model_name: 'all')).models
          descendants = models.to_directed_graph.descendants("project")
          dependent_models = descendants.keys.filter { |k| descendants[k].include?(target_model) }
          ModelsTsv.each_tsv_row(models, dependent_models + [target_model], &block)
        end
      end
    end
  end
end
