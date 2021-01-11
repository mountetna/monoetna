require 'webmock/rspec'
require 'json'
require_relative '../../lib/etna/clients/base_client'

describe 'Base Client class' do
  let(:expired_token) {
    # This is an expired development token and is safe to make public, does not leak anything about production or staging values
    # and cannot be used in a sensitive way.
    # This test depends on breaking down a tok to get expiration info, so it's important we actually set a token
    'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsImZpcnN0IjoiRGV2ZWxvcGVyIiwibGFzdCI6InBlcnNvbiIsInBlcm0iOiJhOnRlc3QtcHJvamVjdCIsImV4cCI6MTAwMH0.cTk-ea5WVpNR3-JXOhC7Z-1n3PL-4NKiZcpha4owr-pZoTtXfM9e3RH5sZW0UTKH7H-DYFdQNV7X1cOLCGTvAOk-3nfOCYHfAfM1sSlyEBEPw3E0dtNCO4APdK_Pz_yooO1yJMQHW0Af4MxTCku4Lu2004JnMsi-hkFATBmNlRprkdaUYHhaSW5rjQ3MSngvku88etZX1-2tFpi_q2FGBmj6_7GfXgxFESmWsGffci4uWceNzntkUIuO_r9xwSsTel99HNvCOnl39YqmnNHx-0uA9BrmyqZGI769f--22tqkpMe_ri-L5pGPZPjcKW2MQhDhl5mSIxCM37SAWsX1Dg'
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