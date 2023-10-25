describe FlagsController do
  include Rack::Test::Methods

  let(:test_flags) { { mythology: "greek", contains_extra_labor: "true" } }

  def app
    OUTER_APP
  end

  def create_flags
    test_flags.each do |k,v|
      create(:flag, project_name: "labors", flag_name: k, value: v)
    end
  end

  context 'get API' do
    it 'shows available flags for authorized users' do
      create_flags
      auth_header(:viewer)
      get('/flags/labors')
      expect(last_response.status).to eq(200)
      expect(json_body[:flags]).to eq(test_flags)
    end

    it 'complains if you are not authorized to view flags' do
      create_flags
      get('/flags/labors')
      expect(last_response.status).to eq(401)
    end
  end

  context 'set API' do

    context 'creates and updates' do

      it 'correctly updates a flag' do
        create_flags
        auth_header(:admin)
        json_post('/flags/labors', flags: [{contains_extra_labor: "false"}])
        get('/flags/labors')
        expect(json_body[:flags][:contains_extra_labor]).to eq("false")
        expect(json_body[:flags][:mythology]).to eq("greek")
      end

      it 'correctly creates flags' do
        create_flags
        auth_header(:admin)
        json_post('/flags/labors', flags: [{is_deprecated: "false"}, {is_rainy: "true"}])
        get('/flags/labors')
        expect(json_body[:flags][:is_deprecated]).to eq("false")
        expect(json_body[:flags][:is_rainy]).to eq("true")
        expect(json_body[:flags][:contains_extra_labor]).to eq("true")
        expect(json_body[:flags][:mythology]).to eq("greek")
      end

      it 'correctly creates a flag if none exist' do
        auth_header(:admin)
        json_post('/flags/labors', flags: [{deprecated: "false"}])
        get('/flags/labors')
        expect(json_body[:flags][:deprecated]).to eq("false")
      end

      it 'complains if you are not authorized to update flags' do
        create_flags
        auth_header(:viewer)
        json_post('/flags/labors', flags: [{deprecated: "false"}])
        expect(last_response.status).to eq(403)
      end
    end

    context 'deletes' do

      it 'correctly deletes a flag' do
        create_flags
        auth_header(:admin)
        json_post('/flags/labors', flags: [{mythology: nil}])
        get('/flags/labors')
        expect(json_body.dig(:flags,:mythology)).to be nil
        expect(json_body[:flags][:contains_extra_labor]).to eq("true")
      end

      it 'gracefully handles when flag does not exist' do
        auth_header(:admin)
        json_post('/flags/labors', flags: [{mythology: nil}])
        get('/flags/labors')
        expect(json_body[:flags]).to be_empty
        expect(last_response.status).to eq(200)
      end
    end

  end
end