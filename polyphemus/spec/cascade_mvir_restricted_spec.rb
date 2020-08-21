describe Polyphemus::CascadeMvirPatientWaiverToRestricted do
  let(:command) { Polyphemus.instance.commands[:cascade_mvir_patient_waiver_to_restricted] }

  let(:all_patients) {
    {
       'Dave': { name: 'Dave', restricted: true, consent: 'Full Study' },
       'Sally': { name: 'Sally', restricted: false, consent: 'Full Study' },
       'Dan': { name: 'Dan', restricted: false, consent: 'Initial Waiver' },
       'Fred': { name: 'Fred', restricted: true, consent: 'Initial Waiver' },
       'Mike': { name: 'Mike', restricted: false, consent: 'No Samples/Data' },
       'Danielle': { name: 'Danielle', restricted: true, consent: 'Full Study' },
    }
  }

  let(:all_pools) { ['pool-a', 'pool-b'] }
  let(:restricted_pools) { ['pool-a'] }

  before(:each) do
    @all_updates = []

    WebMock.disable_net_connect!
    stub_request(:post, 'https://magma.test/retrieve')
        .with(body: hash_including({ project_name: 'mvir1', model_name: 'patient',
                                     attribute_names: 'all', record_names: 'all' }))
        .to_return({ body: {
            'models': { 'patient': { 'documents': all_patients } }
        }.to_json })

    stub_request(:post, 'https://magma.test/update')
        .to_return do |request|

      body = StringIO.new(request.body)
      content_length = body.read.length
      body.rewind

      tempfile = Rack::Multipart::Parser::TEMPFILE_FACTORY
      bufsize = Rack::Multipart::Parser::BUFSIZE
      params = Rack::Utils.default_query_parser

      info = Rack::Multipart::Parser.parse body, content_length, request.headers['Content-Type'], tempfile, bufsize, params

      expect(info.params["project_name"]).to eq("mvir1")
      @all_updates << info.params["revisions"]
      { body: '{}' }
    end

    stub_request(:post, 'https://magma.test/query')
        .with(body: hash_including({ project_name: 'mvir1',
                                     query: [ 'cytof',
                                              [ 'timepoint', 'patient', 'restricted', '::true' ],
                                              '::all', 'cytof_pool', '::identifier' ] }))
        .to_return({ body: {
            'answer': restricted_pools.map {|p| [nil, p] }
        }.to_json })

    stub_request(:post, 'https://magma.test/query')
        .with(body: hash_including({ project_name: 'mvir1',
                                     query: [ 'cytof_pool', '::all', '::identifier' ] }))
        .to_return({ body: {
            'answer': all_pools.map {|p| [nil, p] }
        }.to_json })
  end

  after(:each) do
    WebMock.enable_net_connect!
  end


  it 'updates the restricted states correctly for pools and the patients' do
    command.execute

    expect(@all_updates.inject({}) do |acc, n|
      n.keys.each do |k|
        (acc[k] ||= {}).update(n[k])
      end
      acc
    end).to eq({
      "cytof_pool" => {
        "pool-a" => {"restricted" => "true"},
        "pool-b" => {"restricted" => "false"},
      },
      "patient" => {
          "Mike" => {"restricted" => "true"},
          "Dan" => {"restricted" => "true"},
          "Dave" => {"restricted" => "false"},
          "Danielle" => {"restricted" => "false"}
      }
    })
  end
end
