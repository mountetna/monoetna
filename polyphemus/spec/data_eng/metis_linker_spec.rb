describe MetisLinkerJob do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def labors_config(scripts, bucket_name="pics", autolink=true)
    JSON.parse({
      project_name: "labors",
      config_id: 1,
      version_number: 1,
      config: {
        bucket_name: bucket_name,
        autolink: autolink,
        models: {
          victim: {
            scripts: scripts
          }
        }
      }
    }.to_json)
  end

  let(:models) do
    {
      'victim' => {
       'template' => {
        'identifier' => 'name',
        'parent' => 'monster',
        'attributes' => {
          'monster' => { 'attribute_type' => 'parent' },
          'name' =>  { 'attribute_type' =>  'identifier' },
          'species' =>  { 'attribute_type' =>  'string' },
          'photo_deceased' =>  { 'attribute_type' =>  'file' },
          'target_name' =>  { 'attribute_type' =>  'string' }
        }
       }
      },
      'monster' => {
       'template' => {
        'identifier' => 'name',
        'parent' => 'stub',
        'attributes' => {
          'name' =>  { 'attribute_type' =>  'identifier' },
          'victim' =>  { 'attribute_type' =>  'collection' }
        }
       }
      }
    }
  end

  let(:tail) do
      [
        {"type":"file","id":50,"parent_id":nil,"node_name":"ignore.deceased.png","updated_at":"2023-08-03 22:39:17 +0000","file_hash":"0cc175b9c0f1b6a831c399e269772661","archive_id":nil},
        {"type":"file","id":51,"parent_id":40,"node_name":"LABORS-LION-H2-C1.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"file","id":52,"parent_id":40,"node_name":"NO-ID.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"file","id":53,"parent_id":40,"node_name":"nonmatch_deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"parent","id":40,"parent_id":nil,"node_name":"victims","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil},
        {"type":"parent","id":41,"parent_id":40,"node_name":"archived","updated_at":"2024-02-19 22:39:17 +0000","file_hash":nil,"archive_id":nil},
      ]
  end

  let(:config) {
    config = labors_config([
      {
        type: "file",
        folder_path: "victims",
        file_match: "*.deceased*.png",
        attribute_name: "photo_deceased"
      }
    ])
  }

  let(:rules) do
    {
      monster: '^LABORS-LION$',
      victim: '^LABORS-LION-H\d+-C\d+$'
    }
  end

  let(:runtime_config) {
    {
      "commit" => true
    }
  }
  let(:run_id) { "1234567890" }

  before do
    ENV['TOKEN'] = TEST_TOKEN
    ENV['KUBE_ID'] = run_id
    stub_metis_routes
    stub_request(:post, "https://metis.test/labors/tail/pics").to_return(
      body: tail.to_json,
      headers: { 'Content-Type': "application/json" },
    )
    stub_request(:get, "https://magma.test/gnomon/labors/rules").to_return(
      body: { rules: rules }.to_json,
      headers: { 'Content-Type': "application/json" }
    )
    stub_request(:post, "#{MAGMA_HOST}/retrieve").to_return(
      body: { models: models }.to_json,
      headers: { 'Content-Type': "application/json" }
    )
    stub_magma_update_dry_run

    stub_request(:post, "#{POLYPHEMUS_HOST}/api/workflows/labors/run/update/#{run_id}").to_return(body: "{}")
  end

  context 'linking' do
    it 'successfully links records' do
      stub_request(:post, "https://polyphemus.test/api/workflows/labors/run/previous/1").to_return(
        body: { end_time: Time.now.iso8601 }.to_json,
        headers: { 'Content-Type': "application/json" },
      )

      job = MetisLinkerJob.new(config, runtime_config)

      expect{
        context = job.execute
      }.not_to raise_error
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
      expect(WebMock).to have_requested(:post, /#{POLYPHEMUS_HOST}\/api\/workflows\/labors\/run\/update/)
    end

    it 'successfully links records the first time it is run' do
      stub_request(:post, "https://polyphemus.test/api/workflows/labors/run/previous/1").to_return(
        status: 404,
        body: { error: "No such run for config_id 1 and version_number 1" }.to_json,
        headers: { 'Content-Type': "application/json" },
      )

      job = MetisLinkerJob.new(config, runtime_config)

      expect{
        context = job.execute
      }.not_to raise_error
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
      expect(WebMock).to have_requested(:post, /#{POLYPHEMUS_HOST}\/api\/workflows\/labors\/run\/update/)
    end
  end
end
