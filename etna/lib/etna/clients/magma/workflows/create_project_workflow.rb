# Base workflow for setting up a project by a super user.
# 1) Creates the project in janus
# 2) Adds administrator(s) to Janus
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in Magma.

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
      class CreateProjectWorkflow < Struct.new(:magma_client, :janus_client, :project_name, :project_name_full, keyword_init: true)
        def create_janus_project!
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project_name,
            project_name_full: project_name_full
          ))
        end

        def add_janus_user(email, name, role)
          janus_client.add_user(Etna::Clients::Janus::AddUserRequest.new(
            project_name: project_name,
            email: email,
            name: name,
            role: role,
          ))
        end

        def promote_to_administrator(email)
          janus_client.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
            project_name: project_name,
            # email: user['email'],
            email: email,
            role: 'administrator',
          ))
        end

        def refreshed_token
          janus_client.refresh_token(Etna::Clients::Janus::RefreshTokenRequest.new).token
        end

        def update_magma_client_token!
          self.magma_client = Etna::Clients::Magma.new(
            host: self.magma_client.host,
            token: refreshed_token,
            ignore_ssl: self.magma_client.ignore_ssl)
        end

        def create_magma_project!
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::AddProjectAction.new]))
        end

        def create_magma_project_record!
          magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(
            project_name: project_name,
            revisions: {
                'project' => { project_name => { name: project_name } },
            }))
        end

        def setup_janus_project!
          puts "Creating Janus project."
          create_janus_project!
          puts "Done! Adding you as an administrator on the project."
          add_janus_user(user['email'], "#{user['name']}", 'editor')
          promote_to_administrator(user['email'])
          update_magma_client_token!

          puts "Done with setting up the project in Janus!"
        end

        def setup_magma_project!
          puts "Creating the project in Magma."
          create_magma_project!
          puts "Done! Adding your new project record."
          create_magma_project_record!
        end

        def create!
          setup_janus_project!
          setup_magma_project!

          while true
            puts "Add more users? Y/n"
            break unless STDIN.gets.chomp == 'Y'
            puts "User name?"
            name = STDIN.gets.chomp
            puts "Email?"
            email = STDIN.gets.chomp
            puts "Role? (editor/viewer/administrator)"
            role = STDIN.gets.chomp
            puts "Adding #{name} (#{email}) as a #{role}."
            puts "Confirm? Y/n"
            break unless STDIN.gets.chomp == 'Y'

            if role == 'administrator'
              add_janus_user(email, name, 'editor')
              promote_to_administrator(email)
            else
              add_janus_user(email, name, role)
            end
          end

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
