describe Etna::Application do
  before(:each) do
    class TestApp
      include Etna::Application
    end
  end

  after(:each) do
    Object.send(:remove_const, :TestApp)
  end

  describe 'environment based configuration' do
    before(:each) do
      @orig_ENV = {}.update(ENV)
    end
    after(:each) do
      ENV.clear.update(@orig_ENV)
    end

    it 'loads all kinds of keys successfully, merging them with existing items' do
      ENV['TESTAPP__PRODUCTION__DB__PASSWORD'] = 'password'
      ENV['TESTAPP__PRODUCTION__DB__NO_OVERRIDE'] = 'blah'
      ENV['TESTAPP__PRODUCTION__KEY'] = 'value'
      app = TestApp.instance
      app.configure({:production => { :db => {:host => 'm', :no_override => 'thing'}} })

      expect(app.config(:key, :production)).to eql('value')
      expect(app.config(:db, :production)).to eql({
        :no_override => 'thing',
        :host => 'm',
        :password => 'password'
      })
    end
  end
end
