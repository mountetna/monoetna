# Creates a project from a JSON file.
# This workflow:
# 1) Creates the project in Janus.
# 2) Adds the user as a project administrator to Janus.
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in Magma.
# 5) Iterates through the models and attributes and creates those in Magma.

require 'json'
require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative '../../janus/models'
require_relative './json_models'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class ValidateJsonProjectWorkflow < Struct.new(:filepath, keyword_init: true)
        attr_reader :project

        def initialize(**params)
          super({}.update(params))
          @project = Etna::Clients::Magma::JsonProject.new(filepath: filepath)
        end

        def validate
          raise "Project JSON has errors: ", project.errors unless project.valid?
        end
      end
    end
  end
end
