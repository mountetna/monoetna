require 'digest'


describe 'Integration' do
  
  before do

    # Some dummy data for our dummy request
    @request = {

      'directory'=> '/ipi/melanoma/',
      'expires'=> 600,
      'algorithm'=> 'MD5',
      'timestamp'=> 1473887352,
      'type'=> 'solid',
      'user_email'=> 'jason.cater@ucsf.edu',
      'auth_token'=> 'f976c1a092978d1562e6bb7b16d9a873',
      'file_name'=> 'rand.data.64B' + '-ec65f526379d74e8e4dd60a56f37c868'
    }

    @request_str = '/ipi/melanoma/600MD51473887352solidjason.cater@ucsf.eduf976c1a092978d1562e6bb7b16d9a873rand.data.64B-ec65f526379d74e8e4dd60a56f37c868'
    @request_md5_hash = '20e830f373664673b039239ccbfff5d2'
    @request_sha256_hash = '5262862d77fa71d6aa77fb0d22612e6a0ae331db02e60167c869d88fe5398d14'
  end

  describe Utils do

    it 'sorts request parameters for hashing.' do
      
      ordered_request = Utils.generate_request(@request)
      expect(ordered_request.length).to eq Conf::SIGNATURE_ITEMS.length
    end

    it 'transforms an request from an array to a string.' do

      ordered_request = Utils.generate_request(@request)
      request_str = Utils.stringify_request(ordered_request)
      expect(request_str).to eq @request_str
    end

    it 'hashes a request with an chosen alorithm.' do

      ordered_request = Utils.generate_request(@request)
      hash = Utils.sign_request(ordered_request, @request['algorithm']);
      expect(hash).to eq @request_md5_hash
    end

    it 'hashes a request with MD5.' do

      ordered_request = Utils.generate_request(@request)
      hash = Utils.sign_with_MD5(ordered_request)
      expect(hash).to eq @request_md5_hash
    end

    it 'gashes a request with SHA256.' do

      ordered_request = Utils.generate_request(@request)
      hash = Utils.sign_with_SHA256(ordered_request)
      expect(hash).to eq @request_sha256_hash
    end
  end
end
