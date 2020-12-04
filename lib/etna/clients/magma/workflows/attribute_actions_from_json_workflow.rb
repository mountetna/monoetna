# Updates, renames, or adds attributes to a model, from a JSON file.
# This workflow:
# 1) Validates the attributes JSON file.
# 2) Validates that the model(s) and any link attribute models exist in Magma.
# 3) Executes each of the updates. This will have to be direct Magma client calls,
#     since the model_synch workflow only does Adds and does not currently
#     handle other attribute actions.

require 'json'
require 'ostruct'
require_relative './json_validators'
require_relative './json_converters'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class AttributeActionsFromJsonWorkflow < Struct.new(:magma_client, :project_name, :filepath, keyword_init: true)
        attr_reader :actions
        def initialize(**params)
          super({}.update(params))

          actions_json = JSON.parse(File.read(filepath))

          converter = Etna::Clients::Magma::AttributeActionsConverter.new(actions_json)
          @actions = converter.convert

          validator = Etna::Clients::Magma::AttributeActionsValidator.new(
            actions,
            project_models)
          validator.validate

          raise "Attributes JSON has errors:\n  * #{format_errors(validator.errors)}" unless validator.valid?
        end

        def format_errors(errors)
          errors.map { |e| e.gsub("\n", "\n\t") }.join("\n  * ")
        end

        def project_models
          @project_models ||= magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
            project_name: project_name,
            model_name: 'all')).models
        end

        def execute_actions
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: actions))
        end

        def run!
          puts "Executing the attribute actions against Magma."
          execute_actions
          puts "All complete!"
        end
      end
    end
  end
end
