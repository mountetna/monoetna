require 'webmock/rspec'
require 'json'
require 'base64'
require_relative '../../lib/etna/clients/base_client'

describe 'Base Client class' do
  let(:expired_token) {
    params = {
      email: "user@example.com",
      first: "first",
      last: "last",
      exp: 1000
    }
    "something.#{Base64.strict_encode64(params.to_json)}"
  }

  it 'raises exception if the token is expired' do
    expect {
      client = Etna::Clients::BaseClient.new(token: expired_token, host: 'https://polyphemus.test')
    }.to raise_error(RuntimeError, "Your token is expired.")
  end

  it 'creates correctly with an unexpired token' do
    client = Etna::Clients::BaseClient.new(token: TEST_TOKEN, host: 'https://polyphemus.test')
  end
end