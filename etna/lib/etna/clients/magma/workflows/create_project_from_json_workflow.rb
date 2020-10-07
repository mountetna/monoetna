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
      class CreateProjectFromJsonWorkflow < Struct.new(:magma_client, :janus_client, :filepath, keyword_init: true)
        attr_reader :project, :project_name

        def initialize(**params)
          super({}.update(params))
          @project = Etna::Clients::Magma::JsonProject.new(filepath: filepath)
          raise "Project JSON has errors: ", project.errors unless project.valid?

          @project_name = project.name
        end

        def create_janus_project
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project_name,
            project_name_full: project.name_full
          ))
        end

        def update_janus_permissions
          janus_client.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
            project_name: project_name,
            email: user.email,
            role: 'administrator'
          ))
        end

        def refreshed_token
          janus_client.refresh_token(Etna::Clients::Janus::RefreshTokenRequest.new).token
        end

        def update_magma_client_token
          self.magma_client = Etna::Clients::Magma.new(
            host: self.magma_client.host, token: refreshed_token)
        end

        def create_magma_project
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))
        end

        def add_magma_models
          # Have to add models in a specific order, to ensure that
          #   the parent model exists. So we go from project on down...
          # This assumes that any "link" models exist at a higher level.
          project.models_by_depth.each { |model|
            magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
              project_name: project_name,
              actions: [Etna::Clients::Magma::AddModelAction.new(
                model_name: model.name,
                parent_model_name: model.parent_model_name,
                parent_link_type: model.parent_link_type,
                identifier: model.identifier
              )])
            )
          }
        end

        def add_magma_attributes

        end

        def create!
          create_janus_project
          update_janus_permissions
          update_magma_client_token
          create_magma_project
          add_magma_models
          # add_magma_attributes
        end

        def user
          @user ||= janus_client.whoami(Etna::Clients::Janus::WhoAmIRequest.new).user
        end
      end
    end
  end
end
