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

    def in_temp_file(string)
      Tempfile.open do |file|
        file.write(string)
        file.path
      end
    end

    it 'loads all kinds of keys successfully, merging them with existing items' do
      ENV['ETNA__PRODUCTION__DB__PASSWORD_FILE'] = in_temp_file('- password')
      ENV['ETNA__PRODUCTION__DB__NO_OVERRIDE_FILE'] = in_temp_file('blah')
      ENV['ETNA__PRODUCTION__KEY_FILE'] = in_temp_file('123')
      app = TestApp.instance
      app.configure({:production => { :db => {:host => 'm', :no_override => 'thing'}} })

      expect(app.config(:key, :production)).to eql(123)
      expect(app.config(:db, :production)).to eql({
        :no_override => 'thing',
        :host => 'm',
        :password => ['password']
      })
    end
  end
end
