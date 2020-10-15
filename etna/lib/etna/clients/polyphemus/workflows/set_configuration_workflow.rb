# Given an environment (i.e. test, development, staging),
#   and a polyphemus host, fetches the configuration
#   hash and sets it up in a local config file,
#   ~/.etna.json by default.

# Adds a model to a project, from a JSON file.
# This workflow:
# 1) Validates the model JSON file.
# 2) Validates that the model parent and any link attributes exist in Magma.
# 3) Adds the model to the Magma project.
# 4) Adds any attributes to Magma.

require 'json'
require 'yaml'
require 'ostruct'

module Etna
  module Clients
    class Polyphemus
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class SetConfigurationWorkflow < Struct.new(:host, :environment, :config_file, keyword_init: true)

        def initialize(**params)
          super({config_file: '~/etna.yml'}.update(params))
          prompt_for_token unless ENV['TOKEN']
          @polyphemus_client =
        end

        def fetch_configuration
          polyphemus_client.configuration(Etna::Clients::Polyphemus::ConfigurationRequest.new)
        end

        def update_configuration_file(json_config)
          if !File.exist?(config_file)
            File.open(config_file, 'a+')
          end

          etna_config = YAML.load_file(config_file) || {}

          etna_config[environment] = json_config

          File.open(config_file, 'w') { |f| YAML.dump(etna_config, f) }
        end

        def set_configuration
          json_config = fetch_configuration
          update_configuration_file(json_config)
          puts "Updated ~/etna.yml with your #{environment} configuration."
        end
      end
    end
  end
end
