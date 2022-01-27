describe Etna::EnvironmentVariables do
  class EVTestClass
    include Etna::EnvironmentVariables
  end

  let(:root) { {} }
  let(:env) { {} }
  let(:prefix) { 'PRE' }
  let(:sep) { '__' }
  let(:downcase) { true }
  let(:value_mapper) do
    Proc.new do |path, value|
      [path, value]
    end
  end
  let(:loaded) { EVTestClass.new.load_from_env(prefix, root: root, env: env, downcase: downcase, sep: sep, &value_mapper)}

  describe 'loading environment variables' do
    before(:each) do
      env["#{prefix}__A_B_C"] = "123"
      env["#{prefix}__aValueIsYou"] = "abc"
      env["#{prefix}__some_hash"] = { "a" => 456 }
    end

    it 'works in the happy case' do
      expect(loaded).to eql({"a_b_c"=>"123", "avalueisyou"=>"abc", "some_hash"=>{"a"=>456}})
    end

    describe 'with environment variables lacking the prefix' do
      before(:each) do
        env["some_other_thing"] = "blah"
      end

      it 'does not include that value' do
        expect(loaded).to_not include('some_other_thing')
      end
    end


    describe 'with inner, conflicting value' do
      before(:each) do
        env["#{prefix}__some_hash__b"] = 123
      end

      it 'merges on the parent hash' do
        expect(loaded).to eql({"a_b_c"=>"123", "avalueisyou"=>"abc", "some_hash"=>{"a"=>456, "b" => 123}})
      end
    end
  end
end