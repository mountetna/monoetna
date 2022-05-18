describe Magma::UpdateAttributeAction do
  let(:action) { Magma::RemoveAttributeAction.new("labors", action_params) }
  let(:action_params) do
    {
      action_name: "remove_attribute",
      model_name: "monster",
      attribute_name: "habitat",
    }
  end

  describe '#perform' do
    before do
      @original_attribute = Labors::Monster.attributes[action_params[:attribute_name].to_sym].dup
    end

    after do
      # Rollback in memory changes to the attribute
      new_attr = Magma::Attribute.new(@original_attribute.values)
      new_attr.save
      Labors::Monster.attributes[:habitat] = new_attr
    end

    it 'removes the attribute' do
      expect(action.validate).to be_truthy
      expect(action.perform).to eq(true)
      expect(Labors::Monster.attributes[:habitat]).to be_nil
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
  end
end
