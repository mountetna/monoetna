require 'ostruct'
require_relative '../../json_serializable_struct'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Polyphemus
      class ConfigurationRequest
        def map
          []
        end
      end

      class RedcapJobRequest < Struct.new(:model_names, :redcap_tokens, :commit, :project_name, keyword_init: true)
        include JsonSerializableStruct

        def initialize(**params)
          super({model_names: 'all', commit: false}.update(params))
        end

        def to_json
          {
            job_type: Etna::Clients::Polyphemus::JobType::REDCAP,
            job_params: {
              commit: commit,
              model_names: model_names,
              redcap_tokens: redcap_tokens
            }
          }.to_json
        end
      end

      class ConfigurationResponse
        attr_reader :raw

        def initialize(raw = '')
          @raw = raw
        end

        def environment
          @raw.keys.first
        end

        def environments
          @raw.keys
        end

        def environment_configuration(env = environment)
          EnvironmentConfiguration.new(@raw[env])
        end
      end

      class EnvironmentConfiguration
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def magma
          @raw['magma']
        end

        def metis
          @raw['magma']
        end

        def janus
          @raw['magma']
        end

        def timur
          @raw['magma']
        end

        def polyphemus
          @raw['magma']
        end

        def auth_redirect
          @raw['auth_redirect']
        end
      end

      class JobType < String
        include Enum
        REDCAP = JobType.new("redcap")
      end
    end
  end
end