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

      class ConfigurationResponse
        attr_reader :raw

        def initialize(raw = '')
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
      end
    end
  end
end