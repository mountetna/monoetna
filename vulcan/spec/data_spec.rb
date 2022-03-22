require 'fileutils'

describe DataController do
  include Rack::Test::Methods

  let(:storage) { Vulcan::Storage.new }
  before(:each) do
    clear_store
  end

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'rejects a non-user' do
      auth_header(:non_user)
      store("somehash", "myfile", "abc")
      get("/api/#{PROJECT}/data/somehash/myfile")

      expect(last_response.status).to eq(403)
    end

    it 'fetches the data stored at the given hash and filename' do
      auth_header(:viewer)
      path = store("somehash", "myfile.zip", "abc")
      get("/api/#{PROJECT}/data/somehash/myfile.zip")
      expect(last_response.status).to eq(200)

      expect(last_response['X-Sendfile']).to eq(path)
    end

    it 'returns 404 for a non-existent workflow' do
      auth_header(:viewer)
      get("/api/#{PROJECT}/data/somehash/a")

      expect(last_response.status).to eq(404)

      expect(json_body[:error]).to eq('File a for hash somehash in project labors was not found')
    end
  end
end
