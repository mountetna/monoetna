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
        def plan_synchronization(changeset, project, target_model = 'project')
          apply_matrix_constants(changeset)
          workflow = ModelSynchronizationWorkflow.new(target_client: magma_client,
              source_models: changeset.models, renames: changeset.renames,
              target_project: project, plan_only: true)
          workflow.ensure_model_tree(target_model)
          workflow
        end

        def apply_matrix_constants(changeset)
          changeset.models.model_keys.each do |model_name|
            model = changeset.models.model(model_name)
            attributes = model.template.attributes
            attributes.attribute_keys.each do |attribute_name|
              attribute = attributes.attribute(attribute_name)
              next unless (validation = attribute&.validation)
              next unless (value = validation['value'])

              if validation['type'] == 'Array' && value&.first.start_with?(ModelsCsv::COPY_OPTIONS_SENTINEL)
                digest = value&.first.slice((ModelsCsv::COPY_OPTIONS_SENTINEL.length)..-1)
                attribute.validation = { 'type' => 'Array', 'value' => changeset.matrix_constants[digest] }
              end
            end
          end
        end

        def prepare_changeset_from_csv(csv_io, &err_block)
          line_no = 0
          csv_lines = CSV.parse(csv_io, headers: true, header_converters: :symbol)
          changeset = csv_lines.inject(ModelsCsv::ModelsChangeset.new) do |acc, n|
            line_no += 1
            ModelsCsv.apply_csv_row(acc, n) do |err|
              err_block.call("Error detected on line #{line_no + 1}: #{err}")
              return
            end
          end

          self.class.validate_changeset(changeset, &err_block)

          changeset
        end

        def self.validate_changeset(changeset, &err_block)
          changeset.models.model_keys.each do |model_name|
            validator = AddModelValidator.new(changeset.models, model_name)
            validator.validate
            validator.errors.each(&err_block)
          end

          validator = RenamesValidator.new(changeset.models, changeset.renames)
          validator.validate
          validator.errors.each(&err_block)
        end

        def write_models_template_csv(io, project_name, target_model = 'project')
          models = magma_client.retrieve(RetrievalRequest.new(project_name: project_name, model_name: 'all')).models
          descendants = models.to_directed_graph.descendants(target_model)
          csv = CSV.new(io)
          ModelsCsv.each_csv_row(models, [target_model] + descendants.keys) do |row|
            csv << row
          end
        end
      end
    end
  end
end
