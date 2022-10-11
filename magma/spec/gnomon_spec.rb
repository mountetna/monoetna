
require 'json'

describe GnomonController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  VALID_CONFIG={
    tokens: {
      PROJECT: {
        label: "project",
        values: {
          "The Twelve Labors of Hercules": "The Twelve Labors of Hercules"
        }
      },
      PROJ: {
        label: "project",
        values: {
          "LABORS": "The Twelve Labors of Hercules"
        }
      },
      LABOR: {
        label: "labor",
        values: {
          "The Nemean Lion": "The Nemean Lion",
          "The Lernean Hydra": "The Lernean Hydra"
        }
      },
      LAB: {
        label: "labor",
        values: {
          "LION": "The Nemean Lion",
          "HYDRA": "The Lernean Hydra"
        }
      },
      VILL: {
        label: "Village type",
        values: {
          "V": "Village",
          "H": "Hamlet"
        }
      },
      VICT: {
        label: "Victim type",
        values: {
          "S": "Soldier",
          "C": "Civilian"
        }
      },
      SEP: {
        label: "Separator",
        values: {
          "-": "# Separator"
        }
      }

    },
    synonyms: [
      [ "PROJ", "PROJECT" ],
      [ "LAB", "LABOR" ]
    ],
    rules: {
      project: "PROJECT",
      labor: "LABOR",
      village: "PROJ SEP LAB SEP VILL .n",
      victim: ".village SEP VICT .n"
    }
  }

  it 'complains if there is no grammar' do
    auth_header(:viewer)
    get('/gnomon/labors')

    expect(last_response.status).to eq(422)
  end

  it 'gets the most recent grammar' do
    grammar = create(:grammar, project_name: 'labors', version_number: 1, config: {})
    grammar2 = create(:grammar, project_name: 'labors', version_number: 2, config: VALID_CONFIG)
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

    expect(Magma::Gnomon::Grammar.count).to eq(2)

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
      "property '/text' is invalid: error_type=schema",
      "No separator token defined!"
    ])
    expect(Magma::Gnomon::Grammar.count).to eq(0)
  end

  def create_identifier(id, params={})
    identifier = create(
      :identifier, {
        project_name: 'labors',
        author: "Hera|hera@twelve-labors.org",
        identifier: id,
      }.merge(params)
    )
  end

  it 'decomposes an identifier' do
    Timecop.freeze
    grammar = create(:grammar, project_name: 'labors', version_number: 1, config: VALID_CONFIG)
    identifier = create_identifier("The Twelve Labors of Hercules", rule: 'project', grammar: grammar)
    identifier2 = create_identifier("The Nemean Lion", rule: 'labor', grammar: grammar)
    record = create(:project, name: "The Twelve Labors of Hercules")
    auth_header(:viewer)
    get('/gnomon/labors/decompose/LABORS-LION-H2-C1')

    expect(last_response.status).to eq(200)
    expect(json_body).to eq(
      rules: {
        project: {
          name: "The Twelve Labors of Hercules",
          name_created_at: Time.now.iso8601,
          record_created_at: Time.now.iso8601
        },
        labor: {
          name: "The Nemean Lion",
          name_created_at: Time.now.iso8601,
          record_created_at: nil
        },
        victim: {
          name: "LABORS-LION-H2-C1",
          name_created_at: nil,
          record_created_at: nil
        },
        village: {
          name: "LABORS-LION-H2",
          name_created_at: nil,
          record_created_at: nil
        }
      },
      tokens: [
        ["PROJ", "LABORS"],
        ["SEP", "-"],
        ["LAB", "LION"],
        ["SEP", "-"],
        ["VILL", "H"],
        ["village_counter", "2"],
        ["SEP", "-"],
        ["VICT", "C"],
        ["victim_counter", "1"]
      ]
    )
    Timecop.return
  end

  it 'does not decompose an invalid identifier' do
    grammar = create(:grammar, project_name: 'labors', version_number: 1, config: VALID_CONFIG)
    auth_header(:viewer)
    get('/gnomon/labors/decompose/LABORS-LOON-H2-C1')

    expect(last_response.status).to eq(422)
    expect(json_body[:error]).to eq("Could not decompose identifier LABORS-LOON-H2-C1 for labors")
  end
end
