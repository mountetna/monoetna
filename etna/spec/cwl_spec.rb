describe Etna::Cwl do
  describe 'AnyLoader::ANY' do
    it 'works' do
      val = URI.parse('nonsense')
      expect(Etna::Cwl::AnyLoader::ANY.load(val)).to be(val)
    end
  end

  describe 'SourceLoader' do
    it 'works' do
      expect(Etna::Cwl::SourceLoader.new.load('abc/def/ghi')).to eql(['abc', 'def/ghi'])
      expect(Etna::Cwl::SourceLoader.new.load('abcghi')).to eql([:primary_inputs, 'abcghi'])
      expect(Etna::Cwl::SourceLoader.new.load(['abc', 'def/ghi'])).to eql(['abc', 'def/ghi'])
      expect{Etna::Cwl::SourceLoader.new.load(['abc', 'def', 'ghi'])}.to raise_error
      expect{Etna::Cwl::SourceLoader.new.load(['abc', 2])}.to raise_error
    end
  end

  describe 'TypedDSLLoader::WITH_UNIONS_TYPE_LOADER' do
    def load_with_type(type, val)
      Etna::Cwl::RecordType::Field.type_loader(
          Etna::Cwl::TypedDSLLoader::WITH_UNIONS_TYPE_LOADER.load(type)
      ).load(val)
    end

    it 'can handle record types' do
      record_type = {
          "type" => "record",
          "fields" => {
              "a" => "string",
              "b" => {
                  "type" => {"type" => "enum",
                  "symbols" => ["a", "b", "c"], }
              },
              "c" => {"type" => [{
                  "type" => "record",
                  "fields" => {
                      "g" => "boolean?"
                  },
              }, "null"] }
          }
      }

      value = { 'a' => "hello", "b" => "c", "c" => { "g" => true } }
      expect(load_with_type(record_type, value)).to eql(value)
    end

    it 'can handle explicit union types' do
      expect(load_with_type(['string', 'int'], 145)).to eql(145)
      expect(load_with_type(['string', 'int'], '145')).to eql('145')
      expect{load_with_type(['string', 'int'], 1.45)}.to raise_error
    end

    it 'can handle primitives in various dsl forms' do
      expect(load_with_type('int', 145)).to eql(145)
      expect { load_with_type('int', nil) }.to raise_error
      expect(load_with_type('int?', nil)).to eql(nil)
      expect(load_with_type('int[]', [1, 2, 3])).to eql([1, 2, 3])
      expect{load_with_type('int[]', [1, "a", 3])}.to raise_error
      expect(load_with_type('int[]?', nil)).to be_nil

      expect(load_with_type('string', "abc")).to eql("abc")
      expect { load_with_type('string', 154) }.to raise_error
      expect(load_with_type('boolean', false)).to eql(false)
      expect(load_with_type('boolean', true)).to eql(true)
      expect{load_with_type('boolean', nil)}.to raise_error
      expect{load_with_type('null', 154)}.to raise_error
      expect(load_with_type('null', nil)).to be_nil
      expect(load_with_type('float', 1.45)).to eql(1.45)
      expect(load_with_type('double', 1.45)).to eql(1.45)
      expect(load_with_type('long', 145)).to eql(145)
    end
  end
end
