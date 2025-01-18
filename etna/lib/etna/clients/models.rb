module Etna
  module Clients
    class Response
      attr_reader :raw

      def self.params(*params)
        params.each do |param|
          self.define_method param do
            @raw[param]
          end
        end
      end

      def initialize(raw = {})
        @raw = raw
      end
    end
  end
end
