module Etna::Clients
  describe Magma::AttributesCsv do
    describe '#serialize_attr/#deserialize' do
      Magma::AttributeType.entries.each do |attribute_type|
        next if [
                    Magma::AttributeType::COLLECTION,
                    Magma::AttributeType::CHILD,
                    Magma::AttributeType::TABLE,
                ].include?(attribute_type)

        # Checks both that all attribute types are covered in the case switch, and that they handle the empty string
        # case which is common.
        it "supports deserializing an empty string for #{attribute_type}" do
          attribute = Magma::Attribute.new({
              'attribute_type' => attribute_type,
          })

          deserialized = Magma::AttributesCsv::Importer.deserialize(attribute, '', nil)
          serialized = Magma::AttributesCsv::Exporter.serialize_attr(attribute, deserialized)

          if attribute_type == Magma::AttributeType::BOOLEAN
            expect(serialized).to eql('false')
          else
            expect(serialized).to eql('')
          end
        end
      end
    end
  end
end