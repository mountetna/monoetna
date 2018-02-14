describe Metis::SetUid do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    clear_cookies
  end

  it "sets a metis_uid cookie if none exists" do
    header(
      *Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      )
    )

    get('/')
    expect(last_response.headers["Set-Cookie"]).to match(/#{Metis.instance.config(:metis_uid_name)}/)
    expect(rack_mock_session.cookie_jar[Metis.instance.config(:metis_uid_name)]).not_to be_nil
  end

  it "does not set a cookie if it has been set already" do
    header(
      *Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      )
    )

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=abcdef0123"
    get('/')
    expect(last_response.headers["Set-Cookie"]).to be_nil
    expect(rack_mock_session.cookie_jar[Metis.instance.config(:metis_uid_name)]).not_to be_nil
  end
end
