
require 'json'

describe GnomonController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  VALID_CONFIG={
    tokens: {
      "TOK": {
        label: "Token",
        values: {
          "VAL": "Value"
        }
      },
      "TOKEN": {
        label: "Token",
        values: {
          "VAL": "Value"
        }
      }
    },
    synonyms: [
      [ "TOK", "TOKEN" ]
    ],
    rules: {
      "rule": "TOK"
    }
  }

  it 'complains if there is no grammar' do
    auth_header(:viewer)
    get('/gnomon/labors')

    expect(last_response.status).to eq(422)
  end

  it 'gets the most recent grammar' do
    grammar = create(:grammar, project_name: 'labors', version_number: 1, config: {})
    grammar2 = create(:grammar, project_name: 'labors', version_number: 2, config: JSON.parse(VALID_CONFIG.to_json))
    auth_header(:viewer)
    get('/gnomon/labors')

    expect(last_response.status).to eq(200)
    expect(json_body.to_json).to eq(grammar2.to_hash.to_json)
  end

  it 'sets a new grammar' do
    grammar = create(:grammar, project_name: 'labors', version_number: 1, config: {})

    config = VALID_CONFIG
    auth_header(:admin)
    json_post('/gnomon/labors', config: config)

    expect(last_response.status).to eq(200)

    expect(Magma::Grammar.count).to eq(2)

    expect(json_body.slice(:config, :version_number)).to eq(
      config: config,
      version_number: 2
    )
  end

  it 'does not allow an invalid grammar' do
    config = {
      text: "Some content"
    }
    auth_header(:admin)
    json_post('/gnomon/labors', config: config)

    expect(last_response.status).to eq(422)
    expect(json_body).to match_array(errors: [
      "root is missing required keys: tokens, rules",
      "property '/text' is invalid: error_type=schema"
    ])
    expect(Magma::Grammar.count).to eq(0)
  end
end
