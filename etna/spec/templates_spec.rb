describe Etna::Codegen::Templated do
  class DeferringContainer < Etna::Codegen::Context
    def initialize(template, key)
      super(key)
      @template = template
      @values_watch = template.watch(@values)
    end

    def notify(obj)
      @template.notify(obj)
    end
  end

  class MyTemplate < Etna::Codegen::Context
    include ::Etna::Codegen::Templated
    render_with_erb 'test_template.erb', __FILE__

    def a
      @a ||= DeferringContainer.new(self, {})
    end

    def b
      @b ||= DeferringContainer.new(self, {})
    end

    def maybe

    end

    def items
      watch(@values)
    end
  end
end