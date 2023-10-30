describe FlagsController do
  include Rack::Test::Methods

  let(:test_flags) { { mythology: "greek", contains_extra_labor: "true" } }
  let(:rainy_flag) { { is_rainy: "true" }}
  let(:deprecated_flag) { { is_deprecated: "true" }}

  let(:registered_flags) {
    {
      MYTHOLOGY: {
        name: "mythology",
        greek: "greek"
      },
      EXTRA_LABOR: {
        name: "contains_extra_labor",
        true: "true",
        false: "false"
      },
      RAINY: {
        name: "is_rainy",
        true: "true",
        false: "false"
      },
      DEPRECATED: {
        name: "is_deprecated",
        true: "true"
      }
    }
  }

  def app
    OUTER_APP
  end

  def flag_setup
    register_flag(registered_flags[:MYTHOLOGY])
    register_flag(registered_flags[:EXTRA_LABOR])
    create_flags_in_db(test_flags)
  end

  def unregistered_err_msg(flag_name, flag_value)
    "Flag with name: \"#{flag_name}\", value: \"#{flag_value}\" is not registered in the Magma::Flags module."
  end

  context 'get API' do
    it 'shows available flags for authorized users' do
      flag_setup
      auth_header(:viewer)
      get('/flags/labors')
      expect(last_response.status).to eq(200)
      expect(json_body[:flags]).to eq(test_flags)
    end

    it 'complains if you are not authorized to view flags' do
      flag_setup
      get('/flags/labors')
      expect(last_response.status).to eq(401)
    end
  end

  context 'set API' do

    before do
      unregister_flags
    end

    context 'registered Flags' do

      it 'complains if the flag name is not registered in the Flags module' do
        auth_header(:admin)
        json_post('/flags/labors', flags: [rainy_flag])
        expect(JSON.parse(last_response.body)["error"]).to eq(unregistered_err_msg("is_rainy", "true"))
      end

      it 'complains if the value is not registered in the Flag module ' do
        auth_header(:admin)
        register_flag(registered_flags[:RAINY])
        json_post('/flags/labors', flags: [{is_rainy: "doh!"}])
        expect(JSON.parse(last_response.body)["error"]).to eq(unregistered_err_msg("is_rainy", "doh!"))
      end

      it 'allows flags creates when Flags are registered in the Flag module' do
        auth_header(:admin)
        register_flag(registered_flags[:RAINY])
        json_post('/flags/labors', flags: [rainy_flag])
        get('/flags/labors')
        expect(json_body[:flags][:is_rainy]).to eq("true")
      end
    end

    context 'creates and updates' do

      it 'correctly updates a flag' do
        flag_setup
        auth_header(:admin)
        json_post('/flags/labors', flags: [{contains_extra_labor: "false"}])
        get('/flags/labors')
        expect(json_body[:flags][:contains_extra_labor]).to eq("false")
        expect(json_body[:flags][:mythology]).to eq("greek")
      end

      it 'correctly creates flags' do
        flag_setup
        register_flag(registered_flags[:RAINY])
        register_flag(registered_flags[:DEPRECATED])
        auth_header(:admin)
        json_post('/flags/labors', flags: [deprecated_flag, rainy_flag])
        get('/flags/labors')
        expect(json_body[:flags][:is_deprecated]).to eq("true")
        expect(json_body[:flags][:is_rainy]).to eq("true")
        expect(json_body[:flags][:contains_extra_labor]).to eq("true")
        expect(json_body[:flags][:mythology]).to eq("greek")
      end

      it 'correctly creates a flag if none exist' do
        register_flag(registered_flags[:DEPRECATED])
        auth_header(:admin)
        json_post('/flags/labors', flags: [deprecated_flag])
        get('/flags/labors')
        expect(json_body[:flags][:is_deprecated]).to eq("true")
      end

      it 'complains if you are not authorized to update flags' do
        register_flag(registered_flags[:DEPRECATED])
        auth_header(:viewer)
        json_post('/flags/labors', flags: [deprecated_flag])
        expect(last_response.status).to eq(403)
      end
    end

    context 'deletes' do

      it 'correctly deletes a flag' do
        flag_setup
        auth_header(:admin)
        json_post('/flags/labors', flags: [{mythology: nil}])
        get('/flags/labors')
        expect(json_body.dig(:flags,:mythology)).to be nil
        expect(json_body[:flags][:contains_extra_labor]).to eq("true")
      end

      it 'gracefully handles when flag does not exist' do
        register_flag(registered_flags[:RAINY])
        auth_header(:admin)
        json_post('/flags/labors', flags: [rainy_flag])
        expect(last_response.status).to eq(200)
        json_post('/flags/labors', flags: [{is_rainy: nil}])
        get('/flags/labors')
        expect(json_body[:flags]).to be_empty
        expect(last_response.status).to eq(200)
      end
    end
  end
end
