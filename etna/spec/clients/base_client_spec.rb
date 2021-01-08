require 'webmock/rspec'
require 'json'
require_relative '../../lib/etna/clients/base_client'

describe 'Base Client class' do
  let(:expired_token) {
    # This is an expired development token and is safe to make public, does not leak anything about production or staging values
    # and cannot be used in a sensitive way.
    # This test depends on breaking down a tok to get expiration info, so it's important we actually set a token
    'eyJhbGciOiJSUzI1NiJ9.eyJlbWFpbCI6ImRldmVsb3BlckB1Y3NmLmVkdSIsImZpcnN0IjoiRGV2ZWxvcGVyIiwibGFzdCI6Ikxhc3ROYW1lIiwicGVybSI6IkE6aXBpLG12aXIxLHRlc3QtcHJvamVjdDthOmNvbGUyLGNvbGUzLGNvbGU0LGNvbGVfMSxjb2xlX3g1LHRlc3QxLHRlc3QxMCx0ZXN0MTEsdGVzdDEyLHRlc3QxMyx0ZXN0MTQsdGVzdDE1LHRlc3QxNyx0ZXN0MTgsdGVzdDIzLHRlc3QzLHRlc3Q0LHRlc3Q1LHRlc3Q2LHRlc3Q3LHRlc3Q4LHRlc3Q5O2U6YWRtaW5pc3RyYXRpb24iLCJleHAiOjEwMDB9.maIvYt2Z7H27Rjz7ZBGCVFWrjfxbCrR1nzV2dDu-xZOk8qe3SVORmNUKTWljnRmP4m2KgKyQKFRc2Ninqxo3apLhkzMVDRi4i7ZKKOGCa9OIxwieRbe4XP862Q3CZHKVjQGQ83ofcfkb_8lA-SZwIlPigUKJkXYGPLftbospBEBorut-fTojme16Zub7qLHe12kvzWAsMF7_vfFkwQa97v_U8cDF8zneftwISSQc3BvuqbvwVyqiCF_foBTvxx8-p2OJz7T8SDzighGCOTD_K6pTI6TdJwnS0VcXVs_bwwPvPOtWgu5AsCtk8AvYZMiuZbG9nIEOiFE1bncGvuPVkg'
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