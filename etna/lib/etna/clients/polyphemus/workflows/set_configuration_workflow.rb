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
      class SetConfigurationWorkflow < Struct.new(:polyphemus_client, :config_file, keyword_init: true)
        def fetch_configuration
          polyphemus_client.configuration
        end

        def update_configuration_file(**additional_config)
          if !File.exist?(config_file)
            File.open(config_file, 'a+')
          end

          etna_config = YAML.load_file(config_file) || {}
          config = polyphemus_client.configuration
          env = config.environment.to_sym

          etna_config.delete(env)
          env_config = config.environment_configuration.raw.dup
          env_config.update(additional_config)
          etna_config.update({ env => env_config })

          File.open(config_file, 'w') { |f| YAML.dump(etna_config, f) }
          config
        end
      end
    end
  end
end
