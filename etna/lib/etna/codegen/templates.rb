require 'erb'
require_relative 'contexts'

module Etna
  module Codegen
    class StaticTemplateSegment
      def initialize(txt)
        @txt = txt
      end

      def render
        @txt
      end

      def has_update?
        false
      end
    end

    class MagicBuffer
      def initialize
        @buffer = ''
        @parts = []
      end

      def <<(other)
        if other.is_a? String
          @buffer << other
        else
          flush << other
        end
      end

      def flush
        unless @buffer.empty?
          @parts << StaticTemplateSegment.new(@buffer)
          @buffer = ''
        end
        @parts
      end
    end

    class ErbMagicEscape
      def initialize(v)
        @v = v
      end

      def v
        @v
      end

      def to_s
        v # Shhhh, secret.
      end
    end

    module Templated
      def self.included(mod)
        mod.instance_eval do
          class << self
            def render_with_erb(erb, relative = __FILE__)
              erb = ERB.new(File.open(File.expand_path("../#{erb}", relative)).read, nil, '<>', 'self.buffer')
              erb.def_method(self, '_render(items)', relative)
            end
            
          end
        end
      end

      def buffer=(val)
        @buffer
      end

      def buffer
        @buffer
      end

      def rendered(other_templated)
        ErbMagicEscape.new(other_templated)
      end

      def notify(obj)
        @cached = nil
      end

      def has_update?
        @cached.nil? || child_has_update?
      end

      def child_has_update?
        @parts.any? { |p| p.has_update? }
      end

      def render
        return @cached unless has_update?

        if @cached.nil?
          @parts = []
          @buffer = MagicBuffer.new
          if self.is_a? Context
            _render(@values.values)
          else
            _render([])
          end

          @parts = @buffer.flush
        end

        i = 0
        while child_has_update?
          if (i += 1) > 10
            raise "Update loop failed after #{i} iterations, cyclical dependency is likely."
          end

          @parts.each(&:render)
        end

        @cached = @parts.map(&:render).join('')
      end
    end
  end
end

