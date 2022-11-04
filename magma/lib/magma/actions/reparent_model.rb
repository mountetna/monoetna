class Magma
  class ReparentModelAction < BaseAction
    def perform
      # Swap in the parent attribute to the new parent model.
      # We can then have a project-level migration that detects
      #   the foreign key relation has changed (?)
      #    -- it will still exist in
      #   the project.models, but it will not be linked into
      #   the graph.
      swap_parent_attribute

      @errors.empty?
    end

    def target_models
      if model
        [model]
      else
        []
      end
    end

    def validations
      [
        :validate_model_exists,
        :validate_no_data,
        :validate_parent_model_exists,
        :validate_non_cyclical,
        :validate_different_parent
      ]
    end

    def validate_different_parent
      return unless model && current_parent_model
      return unless current_parent_model.model_name.to_s == @action_params[:parent_model_name]

      @errors << Magma::ActionError.new(
        message: 'parent_model_name is already the parent of model_name.',
        source: @action_params.slice(:model_name, :parent_model_name)
      )
    end

    def validate_no_data
      return unless model
      return unless tree_has_data?

      @errors << Magma::ActionError.new(
        message: 'Cannot reparent a model with data records in its tree (including in child models).',
        source: @action_params.slice(:model_name)
      )
    end

    def validate_model_exists
      return if model

      @errors << Magma::ActionError.new(
        message: 'Model does not exist.',
        source: @action_params.slice(:model_name)
      )
    end

    def validate_parent_model_exists
      return if new_parent_model

      @errors << Magma::ActionError.new(
        message: 'Parent model does not exist.',
        source: @action_params.slice(:parent_model_name)
      )
    end

    def validate_non_cyclical
      # The reparenting cannot result in a cyclical graph!
      return unless model
      return unless new_parent_model

      return unless model_tree.include?(new_parent_model)

      @errors << Magma::ActionError.new(
        message: 'Action would lead to a cyclical graph, since parent_model_name is in the model_name tree.',
        source: @action_params.slice(:model_name, :parent_model_name)
      )
    end

    def model
      @model ||= get_model(@action_params[:model_name])
    rescue NameError => e
      raise e unless e.message =~ /^Could not find.*#{@action_params[:model_name]}$/
      nil
    end

    def new_parent_model
      @new_parent_model ||= get_model(@action_params[:parent_model_name])
    rescue NameError => e
      raise e unless e.message =~ /^Could not find.*#{@action_params[:parent_model_name]}$/
      nil
    end

    def current_parent_model
      @current_parent_model ||= get_model(parent_attribute.link_model_name)
    rescue NameError => e
      raise e unless e.message =~ /^Could not find.*#{parent_attribute.link_model_name}$/
      nil
    end

    private

    def child_type_attributes
      [
        Magma::ChildAttribute,
        Magma::CollectionAttribute,
        Magma::TableAttribute
      ]
    end

    def get_model(model_name)
      Magma.instance.get_model(
        @project_name,
        model_name
      )
    end

    def model_tree
      [model].tap do |result|
        queue = [model]
        seen = []

        current_model = queue.shift
        while current_model
          seen << current_model

          children = child_models_of(current_model)

          # Could get duplicate values if a link is present
          for new_child_model in children
            if !seen.include?(new_child_model)
              queue << new_child_model
              result << new_child_model
            end
          end

          current_model = queue.shift
        end
      end
    end

    def tree_has_data?
      # Check model for records as well as child_models.
      model_tree.any? do |model_to_inspect|
        model_has_data?(model_to_inspect)
      end
    end

    def model_has_data?(model_to_inspect)
      retrieval = Magma::Retrieval.new(
        model_to_inspect,
        'all',
        'identifier',
        collapse_tables: true,
        show_disconnected: true,
        user: user,
        restrict: false # need to check all the data, even if the user cannot see it
      )

      retrieval.count > 0
    end

    def user
      @action_params[:user]
    end

    def child_models_of(model_to_inspect)
      [].tap do |result|
        model_to_inspect.attributes.values.select do |attribute|
          child_type_attributes.include?(attribute.class)
        end.map(&:link_model_name).map do |child_model_name|
          result << get_model(child_model_name)
        end
      end
    end

    def child_model_names
      child_models.map(&:model_name)
    end

    def parent_attribute
      model.attributes.values.select do |attribute|
        attribute.is_a?(Magma::ParentAttribute)
      end.first
    end
  end
end
