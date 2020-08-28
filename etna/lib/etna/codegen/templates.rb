require 'erb'

module Etna
  module Codegen
    module Templated
      def self.included(mod)
        mod.instance_eval do
          attr_reader :cur_context

          class << self
            def render_with_erb(erb, relative = __FILE__)
              erb = ERB.new(File.open(File.expand_path("../#{erb}", relative), eotvar='cur_context').read)
              erb.def_method(self, '_render', relative)
            end
          end

          def render
            @parent_context ||= WatchedContext.new(self, nil) { self._render }
            @parent_context.render
          end

          def watch(obj)
            @segments ||= {}
            @segments[obj.object_id] ||= Watched.new(self, obj)
          end

          def notify(obj)
            @segments ||= {}
            watch = @segments[obj.object_id]
            watch.notify if watch
          end
        end
      end

      class WatchedContext
        attr_reader :children, :parent

        def initialize(templatized, parent, &block)
          @templatized = templatized

          @children = nil # Only initialized during render
          @cached = nil

          @render_block = block
          parent.add_child(self) unless parent.nil?
          @parent = parent
        end

        def add_child(child)
          unless @buffer.empty?
            @buffer = ''
            @children << @buffer
          end

          @children << child
        end

        def <<(str)
          @buffer << str
        end

        def render
          if @children.nil?
            @buffer = ''
            @cached = nil
            @children = []

            prev_context = @templatized.cur_context
            @templatized.cur_context = self
            @render_block.call
            @templatized.cur_context = prev_context
          end

          @cached ||= begin
            result = ''
            @children.each do |c|
              if c.respond_to? :render
                result << c.render
              else
                result << c
              end
            end

            result << @buffer
            result
          end
        end

        def update!(re_evaluate_children)
          return if detached?
          @cached = nil

          if re_evaluate_children && @children
            @children.each(&:detach)
            @children = nil
          end

          render

          parent.update!(false) unless parent.nil? || detached?
        end

        def detach
          @parent = :detached
        end

        def detached?
          parent == :detached
        end
      end

      class Watched
        def initialize(templatized, obj)
          @templatized = templatized
          @obj = obj
          @context_heads = Set.new
          @contexts = []
        end

        def each(&block)
          @contexts << WatchedContext.new(@templatized, @templatized.cur_context) do
            @obj.each(&block)
          end
        end

        def notify
          @contexts.delete_if(&:detached?)
          @contexts.to_a.each { |c| c.update! }
        end
      end
    end
  end
end


# <% mo.each %> do |v|
# <% end %>
# Hi
# <% mo.each %> do |v|
# <% end %>