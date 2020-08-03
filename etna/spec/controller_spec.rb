describe Etna::Controller do
  include Rack::Test::Methods
  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server; end
    end
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
  end

  it "requires parameters to be set" do
    Arachne::Server.get('/test') { require_param(:project_name); success('') }

    @app = setup_app(
      Arachne::Server,
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'e:labors'
    ))

    get('/test')
    expect(last_response.status).to eq(422)

    get('/test?project_name=labors')
    expect(last_response.status).to eq(200)
  end

  context '#logging' do
    before(:each) do
      @log_file = 'error.log'
      Timecop.freeze('2000-01-01')
      # fixes the request_id, which otherwise is random
      Etna::Logger.define_method(:rand) do
        0.23456
      end
    end
    after(:each) do
      ::File.unlink(@log_file)
      Timecop.return
      Etna::Logger.remove_method(:rand)
    end

    it 'reports etna errors in the log file' do
      Arachne::Server.get('/test') { raise Etna::Forbidden, 'You cannot do that.' }

      @app = setup_app(
        Arachne::Server,
        [ Etna::TestAuth ],
        test: { log_file: @log_file, hmac_keys: {etna: 'key'} }
      )

      header(*Etna::TestAuth.token_header(
        email: 'janus@two-faces.org',
        perm: 'e:labors'
      ))

      get('/test')
      expect(last_response.status).to eq(403)

      output = <<EOT
# Logfile created on 2000-01-01 00:00:00 +0000 by logger.rb/61378
ERROR:2000-01-01T00:00:00+00:00 8fzmq8 Exiting with 403, You cannot do that.
EOT
      expect(File.read(@log_file)).to eq(output)
    end

    it 'reports unspecified errors in the log file' do
      Arachne::Server.get('/test') { raise 'Something broke.' }

      @app = setup_app(
        Arachne::Server,
        [ Etna::TestAuth ],
        test: { log_file: @log_file, hmac_keys: {etna: 'key'} },
      )

      header(*Etna::TestAuth.token_header(
        email: 'janus@two-faces.org',
        perm: 'e:labors'
      ))

      get('/test')
      expect(last_response.status).to eq(500)
      output = <<EOT
# Logfile created on 2000-01-01 00:00:00 +0000 by logger.rb/61378
ERROR:2000-01-01T00:00:00+00:00 8fzmq8 Caught unspecified error
ERROR:2000-01-01T00:00:00+00:00 8fzmq8 Something broke.
EOT
      # it reports the error
      log_contents = File.foreach(@log_file).to_a
      expect(log_contents[0..2].join).to eq(output)

      # it reports backtraces
      BACKTRACE = %r!(/[^/]+)+:[0-9]+:in `.*'!
      expect(log_contents[3..-1]).to all( match(BACKTRACE) )
    end
  end
end
