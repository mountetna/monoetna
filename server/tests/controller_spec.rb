require 'rack'
require 'json'
require 'digest'

require './server/tests/dummy_request'
require './server/conf'
require './server/utils'
require './server/controller'

describe 'Integration' do
  
  before do
    
    @controller = Controller.new
    @request = DummyRequest.new
    @request_md5_hash = '20e830f373664673b039239ccbfff5d2'
  end

  describe Controller do

    it 'generates a signature based upon the request.' do

      req = @request.POST
      signature = @controller.generate_signature(req)
      expect(signature).to eq @request_md5_hash
    end

    it 'responds to a valid upload request.' do

      response = @controller.upload(@request)
      response_body = JSON.parse(response.body[0])
      expect(response_body['success']).to eq true
      expect(response_body['signature']).to eq @request_md5_hash
    end

    it 'responds to an invalid upload request.' do
      
      @request.munge # screw up the request data to get an 'invalid' response
      response = @controller.upload(@request)
      response_body = JSON.parse(response.body[0])
      expect(response_body['success']).to eq false
    end
  end
end