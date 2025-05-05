require 'webmock/rspec'

def stub_event_log(project_name='labors')
  polyphemus_host = Etna::Application.instance.config(:polyphemus)[:host]
  WebMock.stub_request(:post, /#{polyphemus_host}\/api\/log\/#{project_name}/).
    to_return(status: 200, body: { log: 1 }.to_json, headers: {'Content-Type': 'application/json'})

  WebMock.stub_request(:options, polyphemus_host).
    to_return({
    status: 200,
    headers: { 'Content-Type': "application/json" },
    body: [
      { :method => "POST", :route => "/api/log/:project_name/write", :name => "log_write", :params => ["project_name"] }
    ].to_json
  })
end
