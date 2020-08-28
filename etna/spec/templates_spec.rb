describe Etna::Codegen::Templated do
  class Template < Etna::Codegen::Context
    include Etna::Codegen::Templated
    render_with_erb 'test_template_self.erb', __FILE__

    def sub_items
      @sub_items ||= TemplateEach.new(a: 1)
    end

    def helper(v)
      "helped #{v}\n"
    end
  end

  class TemplateEach < Etna::Codegen::Context
    include Etna::Codegen::Templated
    render_with_erb 'test_template_each.erb', __FILE__

    def items
      @values.values
    end

    def add_one(v)
      sub_context(key: v) { |key| v }
    end
  end

  it 'works' do
    expect(Template.new(blah: 1).render).to eql("  a
  b
  c
helped t
Hello
World

Thing

")
  end
end