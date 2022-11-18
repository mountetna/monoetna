describe Magma::RemoveLinkAction do
  let(:action) { Magma::RemoveLinkAction.new("labors", action_params) }
  let(:action_params) do
    {
      action_name: "remove_link",
      model_name: "monster",
      attribute_name: "habitat",
    }
  end

  describe '#perform' do
    context 'when column name matches link attribute name' do
      context 'from the link' do
        before do
          @original_attribute = Labors::Monster.attributes[action_params[:attribute_name].to_sym].dup
          @original_reciprocal = Labors::Habitat.attributes[:monster].dup
        end

        after do
          # Rollback in memory changes to the attribute
          Labors::Monster.attributes[action_params[:attribute_name].to_sym] = @original_attribute
          Labors::Habitat.attributes[:monster] = @original_reciprocal
        end

        it 'removes the link' do
          expect(action.validate).to be_truthy
          expect(action.perform).to eq(true)
          expect(Labors::Monster.attributes[action_params[:attribute_name].to_sym]).to be_nil
          expect(Labors::Habitat.attributes[:monster]).to be_nil
        end
      end

      context 'from the reciprocal' do
        let(:action_params) do
          {
            action_name: "remove_link",
            model_name: "habitat",
            attribute_name: "monster",
          }
        end

        before do
          @original_attribute = Labors::Habitat.attributes[action_params[:attribute_name].to_sym].dup
          @original_reciprocal = Labors::Monster.attributes[:habitat].dup
        end

        after do
          # Rollback in memory changes to the attribute
          Labors::Habitat.attributes[action_params[:attribute_name].to_sym] = @original_attribute
          Labors::Monster.attributes[:habitat] = @original_reciprocal
        end

        it 'removes the link' do
          expect(action.validate).to be_truthy
          expect(action.perform).to eq(true)
          expect(Labors::Monster.attributes[:habitat]).to be_nil
          expect(Labors::Habitat.attributes[:monster]).to be_nil
        end
      end
    end

    context 'when column name is different than link attribute name' do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "monster",
          attribute_name: "habitat_name",
        }
      end

      before do
        @original_attribute = Labors::Monster.attributes[:habitat].dup
        @original_reciprocal = Labors::Habitat.attributes[:monster].dup

        rename_action = Magma::RenameAttributeAction.new(
          "labors",
          {
            action: "rename_attribute",
            model_name: "monster",
            attribute_name: "habitat",
            new_attribute_name: "habitat_name"
          }
        )
        rename_action.perform
      end

      after do
        # Rollback in memory changes to the attribute
        Labors::Monster.attributes[:habitat] = @original_attribute
        Labors::Habitat.attributes[:monster] = @original_reciprocal
      end

      it 'removes the attribute' do
        expect(Labors::Monster.attributes[:habitat_name].column_name.to_s).not_to eq("habitat_name")

        expect(action.validate).to be_truthy
        expect(action.perform).to eq(true)
        expect(Labors::Monster.attributes[:habitat_name]).to be_nil
        expect(Labors::Habitat.attributes[:monster]).to be_nil
      end
    end
  end

  describe '#validate' do
    context "when the link attribute name is the identifier" do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "monster",
          attribute_name: "name",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute is not part of a link")
      end
    end

    context "when there's not attribute with attribute_name" do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "monster",
          attribute_name: "not_an_attribute",
        }
      end

      it 'captures an link attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute does not exist")
      end
    end

    context 'when is not a link attribute' do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "monster",
          attribute_name: "species",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute is not part of a link")
      end
    end

    context 'when is a non-link collection attribute' do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "monster",
          attribute_name: "victim",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute is not part of a link")
      end
    end

    context 'when is a parent attribute' do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "victim",
          attribute_name: "monster",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute is not part of a link")
      end
    end

    context 'when is a table attribute' do
      let(:action_params) do
        {
          action_name: "remove_link",
          model_name: "labor",
          attribute_name: "prize",
        }
      end

      it 'captures an attribute error' do
        expect(action.validate).to eq(false)
        expect(action.errors.first[:message]).to eq("Attribute is not part of a link")
      end
    end
  end
end
