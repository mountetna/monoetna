# Bucket permissions are tested using folder#list, but should apply to any
# endpoint calling require_bucket (that is, most of them)

describe Metis::Bucket do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'forbids users from accessing buckets by project role' do

    bucket = create( :bucket, project_name: 'athena', name: 'my_bucket', access: 'editor', owner: 'metis')

    # the admin is always allowed
    token_header(:admin)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(200)

    # the editor is allowed here
    token_header(:editor)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(200)

    # the viewer is forbidden
    token_header(:viewer)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(403)
  end

  it 'forbids users from accessing buckets by email id' do
    bucket = create(:bucket, project_name: 'athena', name: 'my_bucket', access: 'athena@olympus.org', owner: 'metis')

    # the admin is always allowed
    token_header(:admin)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(200)

    # the editor (metis@olympus.org) is forbidden here
    token_header(:editor)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(403)

    # the viewer (athena@olympus.org) is allowed
    token_header(:viewer)
    get('/athena/list/my_bucket/')
    expect(last_response.status).to eq(200)
  end
end

describe BucketController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  context '#list' do
    it 'should return a list of buckets for the current project' do
      bucket2 = create( :bucket, project_name: 'athena', name: 'extra', access: 'viewer', owner: 'metis' )

      token_header(:viewer)
      get('/athena/list')

      expect(last_response.status).to eq(200)

      expect(json_body[:buckets]).to eq(Metis::Bucket.all.map(&:to_hash))
    end
  end
end
