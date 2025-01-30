require_relative '../base_client'

# TODO:  In the near future, I'd like to transition to specifying apis via SWAGGER and generating model stubs from the
# common definitions.  For nowe I've written them out by hand here.
module Etna
  module Clients
    class Gnomon < Etna::Clients::BaseClient
      class ProjectRulesResponse
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def rules
          Rules.new(raw['rules'])
        end
      end

      class Rules
        attr_reader :raw

        def initialize(raw = {})
          @raw = raw
        end

        def [](model_name)
          raw[model_name.to_s]
        end
      end
    end
  end
end
