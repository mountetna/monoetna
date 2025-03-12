require 'webmock/rspec'

def stub_event_log
  WebMock.stub_request(:post, /#{Etna::Application.instance.config(:polyphemus)[:host]}\/api\/log\/labors/).
    to_return(status: 200, body: { log: 1 }.to_json, headers: {'Content-Type': 'application/json'})
end
