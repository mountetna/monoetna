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

              if validation['type'] == 'Array' && value&.first.start_with?(Etna::Clients::Magma::ModelsCsv::COPY_OPTIONS_SENTINEL)
                digest = value&.first.slice((Etna::Clients::Magma::ModelsCsv::COPY_OPTIONS_SENTINEL.length)..-1)
                attribute.validation = { 'type' => 'Array', 'value' => changeset.matrix_constants[digest] }
              end
            end
          end
        end

        def prepare_changeset_from_csv(filename: nil, io: nil, &err_block)
          importer = ModelsCsv::Importer.new
          changeset = importer.prepare_changeset(filename: filename, input_io: io, &err_block)
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

        def write_models_template_csv(project_name, target_model = 'project', filename: nil, io: nil)
          models = magma_client.retrieve(RetrievalRequest.new(project_name: project_name, model_name: 'all')).models
          descendants = models.to_directed_graph.descendants(target_model)
          exporter = ModelsCsv::Exporter.new
          exporter.write_models(models, [target_model] + descendants.keys, filename: filename, output_io: io)
        end
      end
    end
  end
end
