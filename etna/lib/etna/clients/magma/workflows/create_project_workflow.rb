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

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class CreateProjectWorkflow < Struct.new(:magma_client, :janus_client, :project_file, keyword_init: true)
        attr_reader :project

        def initialize(**params)
          super({}.update(params))
          @project = Project.new(project_file)
        end

        def create_janus_project
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project.name,
            project_name_full: project.full_name
          ))
        end

        def create!
        end

        def user
          janus_client.whoami(Etna::Clients::Janus::WhoAmIRequest.new).user
        end
      end

      class Project
        def initialize(project_file)
          @raw = JSON.parse(File.read(project_file), symbolize_names: true)
          raise "Missing required key, \"project_name\"" if !@raw.key?(:project_name)
          raise "Missing required key, \"project_name_full\"" if !@raw.key?(:project_name_full)
        end

        def name
          @raw['project_name']
        end

        def full_name
          @raw['project_name_full']
        end
      end
    end
  end
end
