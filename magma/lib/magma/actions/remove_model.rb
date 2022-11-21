class Magma
  class RemoveModelAction < BaseAction
    def perform
      # Remove the reciprocal attributes from the parent model.
      # We can then have a project-level migration that detects
      #   when a model is detached -- it will still exist in
      #   the project.models, but it will not be linked into
      #   the graph.
      remove_reciprocal_attributes

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
        :validate_leaf,
        :validate_no_links
      ]
    end

    def validate_leaf
      return unless model
      return unless model.attributes.values.any? do |attribute|
        child_type_attributes.include?(attribute.class)
      end

      @errors << Magma::ActionError.new(
        message: 'Model is not a leaf model.',
        source: @action_params.slice(:model_name)
      )
    end

    def validate_no_links
      link_type_attributes = [
        Magma::LinkAttribute
      ]

      return unless model
      return unless model.attributes.values.any? do |attribute|
        link_type_attributes.include?(attribute.class)
      end

      # But, links that refer to the same model are okay, i.e.
      #   reference_patient
      link_attributes = model.attributes.values.select do |attribute|
        link_type_attributes.include?(attribute.class)
      end

      return unless link_attributes.any? do |attribute|
        attribute.link_model_name != @action_params[:model_name]
      end

      @errors << Magma::ActionError.new(
        message: 'Model has a link to another model. Use "remove_link" on that attribute, first.',
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

    def model
      @model ||= Magma.instance.get_model(
        @project_name,
        @action_params[:model_name]
      )
    rescue NameError => e
      raise e unless e.message =~ /^Could not find.*#{@action_params[:model_name]}$/
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

    def remove_reciprocal_attributes
      reciprocal_models.each do |reciprocal_attribute, reciprocal_model|
        reciprocal_model.attributes.delete(reciprocal_attribute.attribute_name.to_sym)
        reciprocal_attribute.delete
      end
    end

    def reciprocal_models
      foreign_key_attributes.map do |attribute|
        reciprocal_model = Magma.instance.get_model(
          @project_name,
          attribute.link_model_name
        )

        [
          reciprocal_model.attributes[attribute.link_attribute_name.to_sym],
          reciprocal_model
        ]
      end
    end

    def foreign_key_attributes
      model.attributes.values.select do |attribute|
        attribute.is_a?(Magma::ForeignKeyAttribute)
      end
    end
  end
end
