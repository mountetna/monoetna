# Base workflow for setting up a project by a super user.
# 1) Creates the project in janus
# 2) Adds administrator(s) to Janus
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in .

require 'base64'
require 'json'
require 'ostruct'
require_relative '../models'

module Etna
  module Clients
    class Janus
      class GenerateTokenWorkflow < Struct.new(:janus_client, :email, :project_name, :token_type, :private_key_file, keyword_init: true)
        def generate!
          nonce = janus_client.get_nonce

          unless email
            puts "Email address for #{janus_client.host} account?"
            email = STDIN.gets.chomp
          end

          if use_nonce?
            until private_key_file
              puts "Location of private key file?"
              private_key_file = ::File.expand_path(STDIN.gets.chomp)
              unless File.exists?(private_key_file)
                puts "No such file."
                private_key_file = nil
              end
            end
          end

          if needs_project_name?
            puts "Project name?"
            project_name = STDIN.gets.chomp
          end
          
          token = janus_client.generate_token(token_type, signed_nonce: use_nonce? ? signed_nonce(nonce) : nil, project_name: project_name)

          puts token
        end

        private

        def signed_nonce(nonce)
          private_key = OpenSSL::PKey::RSA.new(File.read(private_key_file))

          txt_to_sign = "#{nonce}.#{Base64.strict_encode64(email)}"

          sig = Base64.strict_encode64(
            private_key.sign(OpenSSL::Digest::SHA256.new,txt_to_sign)
          )

          "#{txt_to_sign}.#{sig}"
        end

        def use_nonce?
          !task_token? || !janus_client.token
        end

        def needs_project_name?
          task_token? && !project_name
        end

        def task_token?
          token_type == 'task'
        end

        def user
          @user ||= JSON.parse(Base64.urlsafe_decode64(magma_client.token.split('.')[1]))
        end
      end
    end
  end
end
