require 'rack'
require 'json'
require 'digest'

require './server/conf'
require './server/utils'
require './server/controller'

describe 'Integration' do
  
  before do

    @controller = Controller.new

    # Some dummy data for our dummy request
    @request = {

      params: {

        directory: '/ipi/melanoma/',
        expires: 600,
        algorithm: 'MD5',
        timestamp: 1473887352,
        type: 'solid',
        user_email: 'jason.cater@ucsf.edu',
        auth_token: 'f976c1a092978d1562e6bb7b16d9a873',
        file_name: 'rand.data.64B' + '-ec65f526379d74e8e4dd60a56f37c868',
        signature: '20e830f373664673b039239ccbfff5d2'
      }
    }

    @request_md5_hash = '20e830f373664673b039239ccbfff5d2'
  end

  describe Controller do

    it 'Generates a signature based upon the request.' do

      signature = @controller.generate_signature(@request)
      expect(signature).to eq @request_md5_hash
    end

    it 'Responds to a valid upload request.' do
    
      response = @controller.upload(@request)
      response_body = JSON.parse(response.body[0])
      expect(response_body['success']).to eq true
      expect(response_body['signature']).to eq @request_md5_hash
    end

    it 'Responds to an invalid upload request.' do
      
      @request[:params][:signature] = ''
      response = @controller.upload(@request)
      response_body = JSON.parse(response.body[0])
      expect(response_body['success']).to eq false
    end
  end
end