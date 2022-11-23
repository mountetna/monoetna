describe EtnaApp::Config::Set do
  let(:etna_config_file) { Tempfile.new }

  before(:each) do
    etna_config_file = self.etna_config_file
    EtnaApp.define_singleton_method(:config_file_path) { etna_config_file.path }
  end

  after(:each) do
    etna_config_file.close!
  end

  include WithEtnaClients
  include WithLogger

  # Set TOKEN before re-recording test
  it 'collects server configuration into a local config file' do
    VCR.use_cassette('set_configuration_spec.e2e') do
      ENV['TOKEN'] ||= TEST_TOKEN
      cmd, args, kwds = EtnaApp.instance.find_command('config', 'set', 'https://polyphemus.development.local', '--ignore-ssl')
      cmd.execute(*args, **kwds)
    end

    contents = YAML.load(File.open(etna_config_file.path) { |f| f.read })
    expect(contents).to eql({
        development: {
            auth_redirect: "https://janus.development.local",
            janus: { host: "https://janus.development.local" },
            magma: { host: "https://magma.development.local" },
            metis: { host: "https://metis.development.local" },
            polyphemus: { host: "https://polyphemus.development.local" },
            timur: { host: "https://timur.development.local" },
            ignore_ssl: true,
            docker: nil
        },
    })

    EtnaApp.instance.configure(contents)
    expect(magma_client).to_not be_nil
  end
end

