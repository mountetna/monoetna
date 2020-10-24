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
      class CreateProjectFromJsonWorkflow < Struct.new(:magma_client, :janus_client, :filepath, keyword_init: true)
        attr_reader :project, :converter

        def initialize(**params)
          super

          user_json = JSON.parse(File.read(filepath))
          magma_json = Etna::Clients::Magma::ProjectConverter.convert_project_user_json_to_magma_json(user_json)
          @project = Etna::Clients::Magma::Project.new(magma_json)

          @validator = Etna::Clients::Magma::ProjectValidator.new(project)
          @validator.validate

          raise "Project JSON has errors: #{@validator.errors}" unless @validator.valid?

          @converter = ProjectConverter.new(project)
          @converter.convert!
        end

        def project_name
          project.raw['project_name']
        end

        def project_name_full
          project.raw['project_name_full']
        end

        def create_janus_project
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project_name,
            project_name_full: project_name_full
          ))
        end

        def add_janus_user
          janus_client.add_user(Etna::Clients::Janus::AddUserRequest.new(
            project_name: project_name,
            email: user['email'],
            name: "#{user['first']} #{user['last']}",
            role: 'editor'
          ))
        end

        def update_janus_permissions
          janus_client.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
            project_name: project_name,
            email: user['email'],
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
          converter.model_tree.each do |model|
            next if model.name == 'project'

            magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
              project_name: project_name,
              actions: [Etna::Clients::Magma::AddModelAction.new(
                model_name: model.name,
                parent_model_name: model.raw['parent_model_name'],
                parent_link_type: model.raw['parent_link_type'],
                identifier: model.template.identifier
              )]))
          end
        end

        def ensure_magma_tree
          # Call the model_synchronization_workflow.ensure_model_tree
          #   on each model to set their attributes in the target Magma.
          workflow = Etna::Clients::Magma::ModelSynchronizationWorkflow.new(
            target_project: project_name,
            target_client: magma_client,
            source_models: project.models
          )
          project.models.model_keys.each do |model_name|
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
          project.models.model_keys.each do |model_name|
            project.models.model(model_name).template.attributes.all.select do |attribute|
              attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
            end.each do |attribute|
              magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
                project_name: project_name,
                actions: [Etna::Clients::Magma::UpdateAttributeAction.new(
                  model_name: model_name,
                  attribute_name: attribute.attribute_name,
                  display_name: attribute.display_name,
                  description: attribute.desc,
                  validation: attribute.validation
                )]))
            end
          end
        end

        def setup_janus_project!
          puts "Creating Janus project."
          create_janus_project
          puts "Done! Adding you as an administrator on the project."
          add_janus_user
          update_janus_permissions
          update_magma_client_token
          puts "Done with setting up the project in Janus!"
        end

        def setup_magma_project!
          puts "Creating the project in Magma."
          create_magma_project
          puts "Done! Creating all the models and attributes in Magma."
          create_magma_models
          ensure_magma_tree
          puts "Done! Adding your new project record."
          create_magma_project_record
          update_magma_attributes
        end

        def create!
          setup_janus_project!
          setup_magma_project!
          puts "All complete!"
          puts "You need to visit Janus to refresh your token."
          puts "You can now log into any app to manage your data."
        end

        def user
          @user ||= JSON.parse(Base64.urlsafe_decode64(magma_client.token.split('.')[1]))
        end
      end
    end
  end
end
