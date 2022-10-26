describe Magma::UpdateAttributeAction do
  let(:action) { Magma::RemoveAttributeAction.new("labors", action_params) }
  let(:action_params) do
    {
      action_name: "remove_attribute",
      model_name: "monster",
      attribute_name: "species",
    }
  end

  describe '#perform' do
    context 'when column name matches attribute name' do
      before do
        @original_attribute = Labors::Monster.attributes[action_params[:attribute_name].to_sym].dup
      end

      after do
        # Rollback in memory changes to the attribute
        Labors::Monster.attributes[action_params[:attribute_name].to_sym] = @original_attribute
      end

      it 'removes the attribute' do
        expect(action.validate).to be_truthy
        expect(action.perform).to eq(true)
        expect(Labors::Monster.attributes[action_params[:attribute_name].to_sym]).to be_nil
      end
    end

    context 'when column name is different than attribute name' do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "monster",
          attribute_name: "species_name",
        }
      end

      before do
        @original_attribute = Labors::Monster.attributes[:species].dup

        rename_action = Magma::RenameAttributeAction.new(
          "labors",
          {
            action: "rename_attribute",
            model_name: "monster",
            attribute_name: "species",
            new_attribute_name: "species_name"
          }
        )
        rename_action.perform
      end

      after do
        # Rollback in memory changes to the attribute
        Labors::Monster.attributes[:species] = @original_attribute
      end

      it 'removes the attribute' do
        expect(Labors::Monster.attributes[:species_name].column_name.to_s).not_to eq("species_name")

        expect(action.validate).to be_truthy
        expect(action.perform).to eq(true)
        expect(Labors::Monster.attributes[:species_name]).to be_nil
      end
    end
  end

  describe '#validate' do
    context "when the attribute name is the identifier" do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "monster",
          attribute_name: "name",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Cannot remove identity attribute")
      end
    end

    context "when there's not attribute with attribute_name" do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "monster",
          attribute_name: "not_an_attribute",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute does not exist")
      end
    end

    context 'when is a link attribute' do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "monster",
          attribute_name: "habitat",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Use remove_link action to remove a link attribute")
      end
    end

    context 'when is a collection attribute' do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "monster",
          attribute_name: "victim",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Use reparent_model action to change this attribute")
      end
    end

    context 'when is a parent attribute' do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "victim",
          attribute_name: "monster",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Use reparent_model action to change this attribute")
      end
    end

    context 'when is a table attribute' do
      let(:action_params) do
        {
          action_name: "remove_attribute",
          model_name: "labor",
          attribute_name: "prize",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Use reparent_model action to change this attribute")
      end
    end
  end
end
