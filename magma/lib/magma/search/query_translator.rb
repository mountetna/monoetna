require_relative '../audit/project_to_template_mapping'
require_relative '../audit/model_hierarchy'

class Magma
  module Search
    class QueryTranslator
      class NotApplicable < StandardError
      end

      # Filter verbs whose argument may name a template model. The argument is
      # replaced with the link attribute from the current model to that
      # model's child model.
      LINK_VERBS = ['::has', '::lacks'].freeze

      def initialize(project_name)
        @project_name = project_name
        @mapping = Magma::ProjectToTemplateMapping.new(project_name)
        @hierarchy = Magma::ModelHierarchy.new(project_name)
      end

      def translate(query)
        check_project_loaded!
        translated = translate_node(query, nil)
        to_plain_json(translated)
      end

      private

      def check_project_loaded!
        raise NotApplicable, "#{@project_name} is not loaded" unless Magma.instance.get_project(@project_name)
      end

      # Magma::Question expects plain JSON (string keys, no symbols).
      def to_plain_json(query)
        JSON.parse(JSON.generate(query))
      end

      # A query is a nested array; each level starts with a model, a verb, or neither.
      def translate_node(node, parent_template)
        return node unless node.is_a?(Array)

        head, *args = node
        if verb?(head)
          translate_verb(head, args, parent_template)
        elsif template_model?(head)
          translate_model(head, args, parent_template)
        else
          translate_args(node, parent_template)
        end
      end

      def translate_args(args, parent_template)
        args.map { |part| translate_node(part, parent_template) }
      end

      def translate_model(template_model, args, parent_template)
        translated = [child_model!(template_model), *translate_args(args, template_model)]
        add_intermediate_models(translated, template_model, parent_template)
      end

      def translate_verb(verb, args, parent_template)
        if link_verb?(verb, args)
          translate_link_verb(verb, args, parent_template)
        else
          [verb, *translate_args(args, parent_template)]
        end
      end

      def translate_link_verb(verb, args, parent_template)
        raise NotApplicable, 'link verb outside a model' unless parent_template

        target_child = child_model!(args.first)
        context_child = child_model!(parent_template)
        link = @hierarchy.link_between(context_child, target_child)
        raise NotApplicable, "no link from #{context_child} to #{target_child}" unless link

        [verb, link[:attribute_name], *translate_args(args.drop(1), parent_template)]
      end

      # parent_path walks from the child model up to the template parent's
      # child model. Steps after the first are models that sit between the
      # two template hops, nearest the child first. Wrap those as extra
      # ::any levels, outermost nearest the parent.
      def add_intermediate_models(translated, template_model, parent_template)
        return translated if parent_template.nil?

        child_name = child_model!(template_model)
        parent_child = child_model!(parent_template)
        path = @hierarchy.parent_path(child_name, parent_child)
        raise NotApplicable, "no parent path from #{child_name} to #{parent_child}" if path.nil?

        path.drop(1).map { |step| step[:model] }.reverse.each do |extra|
          translated = [extra, translated, '::any']
        end
        translated
      end

      def verb?(part)
        part.is_a?(String) && part.start_with?('::')
      end

      def link_verb?(verb, args)
        LINK_VERBS.include?(verb) && template_model?(args.first)
      end

      def template_model?(name)
        name.is_a?(String) && @mapping.template_model_names.include?(name)
      end

      def child_model!(template_model)
        @mapping.child_model_name(template_model) || raise(NotApplicable, "#{template_model} is not mapped")
      end
    end
  end
end
