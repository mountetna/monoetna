# Creates a project from a JSON file.
# This workflow:
# 1) Creates the project in Janus.
# 2) Adds the user as a project administrator to Janus.
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in Magma.
# 5) Iterates through the models and attributes and creates those in Magma.

require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative '../../janus/models'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class CreateProjectWorkflow < Struct.new(:magma_client, :janus_client, :project_file, keyword_init: true)
        def create_janus_project
          janus_client.add_project()
        end

        def create!
        end
      end
    end
  end
end
