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
          raise "Project JSON has errors: #{project.errors}" unless project.valid?

          @project_name = project.name
        end

        def create_janus_project
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project_name,
            project_name_full: project.name_full
          ))
        end

        def add_janus_user
          janus_client.add_user(Etna::Clients::Janus::AddUserRequest.new(
            project_name: project_name,
            email: user.email,
            name: "#{user.first} #{user.last}",
            role: 'editor'
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

        def create_magma_models
          # Technically ensure_magma_tree should create models as needed,
          #   but things seem to go more smoothly if we manually
          #   create the models before we try to ensure the tree.
          project.model_tree.each do |model|
            magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
              project_name: project_name,
              actions: [Etna::Clients::Magma::AddModelAction.new(
                model_name: model.name,
                parent_model_name: model.parent_model_name,
                parent_link_type: model.parent_link_type,
                identifier: model.identifier
              )]))
          end
        end

        def magma_models
          # Convert each of our JSON models into an instance of
          #   Magma Model + Attributes
          @magma_models ||= project.get_magma_models
        end

        def ensure_magma_tree
          # Call the model_synchronization_workflow.ensure_model_tree
          #   on each model to set their attributes in the target Magma.
          workflow = Etna::Clients::Magma::ModelSynchronizationWorkflow.new(
            target_project: project.name,
            target_client: magma_client,
            source_models: magma_models
          )
          magma_models.model_keys.each do |model_name|
            workflow.ensure_model_tree(model_name)
          end
        end

        def create_magma_project_record
          revision = {}
          revision[project_name] = {
            name: project_name
          }
          magma_client.update(Etna::Clients::Magma::UpdateRequest.new(
            project_name: project_name,
            revisions: {
              project: revision
            }))
        end

        def update_magma_attributes
          # Sometimes an identifier attribute might have
          #   validations or other settings included in the
          #   model definition. So we have to "update" those
          #   instead of adding those (which occurs in
          #   ensure_magma_tree).
          magma_models.model_keys.each do |model_name|
            magma_models.model(model_name).template.attributes.all.select do |attribute|
              attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
            end.each do |attribute|
              magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
                project_name: project_name,
                actions: [Etna::Clients::Magma::UpdateAttributeAction.new(
                  model_name: model_name,
                  attribute_name: attribute.name,
                  display_name: attribute.display_name,
                  description: attribute.desc,
                  validation: attribute.validation
                )]))
            end
          end
        end

        def create!
          create_janus_project
          add_janus_user
          update_janus_permissions
          update_magma_client_token
          create_magma_project
          create_magma_models
          ensure_magma_tree
          create_magma_project_record
          update_magma_attributes
        end

        def user
          @user ||= janus_client.whoami(Etna::Clients::Janus::WhoAmIRequest.new).user
        end
      end
    end
  end
end
