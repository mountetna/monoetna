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
      class AddProjectModelsWorkflow < Struct.new(:magma_client, keyword_init: true)
        def synchronize_to_server(models, project, target_model = 'project', &update_block)
          workflow = ModelSynchronizationWorkflow.new(target_client: magma_client, source_models: models, target_project: project, update_block: update_block)
          workflow.ensure_model_tree(target_model)
        end

        def prepare_models_from_csv(csv_io, &err_block)
          line_no = 0
          csv_lines = CSV.parse(csv_io, headers: true, header_converters: :symbol)
          models = csv_lines.inject(Models.new) do |acc, n|
            line_no += 1
            ModelsCsv.apply_csv_row(acc, n) do |err|
              err_block.call("Error detected on line #{line_no + 1}: #{err}")
              return
            end
          end

          models.model_keys.each do |model_name|
            validator = AddModelValidator.new(models, model_name)
            validator.validate
            validator.errors.each(&err_block)
          end

          models
        end

        def write_models_templats_csv(csv, project_name, target_model = 'project')
          models = magma_client.retrieve(RetrievalRequest.new(project_name: project_name, model_name: 'all')).models
          descendants = models.to_directed_graph.descendants(target_model)
          csv = CSV.new(csv)
          ModelsCsv.each_csv_row(models, [target_model] + descendants.keys) do |row|
            csv << row
          end
        end
      end
    end
  end
end
