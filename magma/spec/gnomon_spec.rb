
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

  context 'counter API' do
    it 'throws exception with no grammar defined' do
      expect(Magma::Gnomon::Identifier.count).to eq(0)
      auth_header(:admin)
      post('/gnomon/labors/increment/victim/LABORS-LION-H2-C')
      expect(last_response.status).to eq(422)
      expect(Magma::Gnomon::Identifier.count).to eq(0)
    end

    context 'with valid grammar' do
      before(:each) do
        @grammar = create(:grammar, project_name: 'labors', version_number: 1, config: VALID_CONFIG)
      end

      context 'generates the next identifier when' do
        it 'none exist' do
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-C')
          expect(last_response.status).to eq(200)
          expect(last_response.body).to eq("1")
        end

        it 'does not create the identifier' do
          expect(Magma::Gnomon::Identifier.count).to eq(0)
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-C')
          expect(last_response.status).to eq(200)
          expect(last_response.body).to eq("1")
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'sequence exists' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-C')
          expect(last_response.status).to eq(200)
          expect(last_response.body).to eq("2")
        end

        it 'sequence in other token value exists' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-S')
          expect(last_response.status).to eq(200)
          expect(last_response.body).to eq("1")
        end

        it 'change in parent token counter' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H1-C')
          expect(last_response.status).to eq(200)
          expect(last_response.body).to eq("1")
        end
      end

      context 'throws exception when' do
        it 'user is not an administrator' do
          auth_header(:viewer)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-C')
          expect(last_response.status).to eq(403)
        end

        it 'rule does not exist' do
          auth_header(:admin)
          post('/gnomon/labors/increment/habitat/LABORS-MARSH')
          expect(last_response.status).to eq(422)
        end

        it 'identifier_root does not match rule' do
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-Q')
          expect(last_response.status).to eq(422)
          post('/gnomon/labors/increment/victim/LABORS-PARROT-H')
          expect(last_response.status).to eq(422)
        end

        it 'but does not create identifier' do
          expect(Magma::Gnomon::Identifier.count).to eq(0)
          auth_header(:admin)
          post('/gnomon/labors/increment/victim/LABORS-LION-H2-Q')
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end
      end
    end
  end

  context 'list API' do

    it 'throws exception when no grammar for project' do
      auth_header(:viewer)
      get('/gnomon/labors/list/victim')
      expect(last_response.status).to eq(422)
    end

    context 'with grammar' do
      before(:each) do
        @grammar = create(:grammar, project_name: 'labors', version_number: 1, config: VALID_CONFIG)
      end

      it 'returns all identifiers with no regex param' do
        identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
        auth_header(:viewer)
        get('/gnomon/labors/list/victim')
        expect(last_response.status).to eq(200)

        expect(json_body.length).to eq(1)
        expect(json_body.first).to include({
          identifier: "LABORS-LION-H2-C1",
          author: "Hera|hera@twelve-labors.org"})
        expect(json_body.first[:name_created_at]).not_to eq(nil)
        expect(json_body.first[:record_created_at]).to eq(nil)
      end

      it 'correctly filters using the regex param' do
        identifier1 = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
        identifier2 = create_identifier("LABORS-LION-H2-C2", rule: 'victim', grammar: @grammar)

        auth_header(:viewer)
        get('/gnomon/labors/list/victim?regex=C2')
        expect(last_response.status).to eq(200)

        expect(json_body.length).to eq(1)
        expect(json_body.first).to include({
          identifier: "LABORS-LION-H2-C2",
          author: "Hera|hera@twelve-labors.org"})

        get('/gnomon/labors/list/victim?regex=LION')
        expect(last_response.status).to eq(200)

        expect(json_body.length).to eq(2)
        expect(json_body.map { |id| id[:identifier] }).to match_array([
          "LABORS-LION-H2-C1",
          "LABORS-LION-H2-C2"
        ])
      end

      context 'restricted records' do
        it 'non-privileged user does not see record_created_at times' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)

          project = create(:project, name: 'The Twelve Labors of Hercules')
          lion = create(:labor, :lion, project: project)

          lion_monster = create(:monster, :lion, labor: lion)

          victim = create(:victim, name: 'LABORS-LION-H2-C1', monster: lion_monster, country: 'Italy', restricted: true)

          auth_header(:viewer)
          get('/gnomon/labors/list/victim')
          expect(last_response.status).to eq(200)

          expect(json_body.length).to eq(1)
          expect(json_body.first).to include({
            identifier: "LABORS-LION-H2-C1",
            author: "Hera|hera@twelve-labors.org"})
          expect(json_body.first[:name_created_at]).not_to eq(nil)
          expect(json_body.first[:record_created_at]).to eq(nil)
        end

        it 'privileged user does see record_created_at times' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)

          project = create(:project, name: 'The Twelve Labors of Hercules')
          lion = create(:labor, :lion, project: project)

          lion_monster = create(:monster, :lion, labor: lion)

          victim = create(:victim, name: 'LABORS-LION-H2-C1', monster: lion_monster, country: 'Italy', restricted: true)

          auth_header(:privileged_editor)
          get('/gnomon/labors/list/victim')
          expect(last_response.status).to eq(200)

          expect(json_body.length).to eq(1)
          expect(json_body.first).to include({
            identifier: "LABORS-LION-H2-C1",
            author: "Hera|hera@twelve-labors.org"})
          expect(json_body.first[:name_created_at]).not_to eq(nil)
          expect(json_body.first[:record_created_at]).not_to eq(nil)
        end
      end

      it 'correctly supplies record creation times from Magma for unrestricted records' do
        identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)

        project = create(:project, name: 'The Twelve Labors of Hercules')
        lion = create(:labor, :lion, project: project)

        lion_monster = create(:monster, :lion, labor: lion)

        victim = create(:victim, name: 'LABORS-LION-H2-C1', monster: lion_monster, country: 'Italy')

        auth_header(:viewer)
        get('/gnomon/labors/list/victim')
        expect(last_response.status).to eq(200)

        expect(json_body.length).to eq(1)
        expect(json_body.first).to include({
          identifier: "LABORS-LION-H2-C1",
          author: "Hera|hera@twelve-labors.org"})
        expect(json_body.first[:name_created_at]).not_to eq(nil)
        expect(json_body.first[:record_created_at]).not_to eq(nil)
      end

      context 'throws exception when' do
        it 'invalid rule name provided' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
          auth_header(:viewer)
          get('/gnomon/labors/list/alias')
          expect(last_response.status).to eq(422)
        end
      end
    end
  end

  context 'generate API' do

    it 'throws exception when no grammar for project' do
      auth_header(:admin)
      post('/gnomon/labors/generate/victim/LABORS-LION-H2-C1')
      expect(last_response.status).to eq(422)
      expect(Magma::Gnomon::Identifier.count).to eq(0)
    end

    context 'with grammar' do
      before(:each) do
        @grammar = create(:grammar, project_name: 'labors', version_number: 1, config: VALID_CONFIG)
      end

      it 'creates the identifier' do
        identifier = "LABORS-LION-H2-C1"
        expect(Magma::Gnomon::Identifier.count).to eq(0)
        auth_header(:admin)
        post("/gnomon/labors/generate/victim/#{identifier}")
        expect(last_response.status).to eq(200)
        expect(json_body[:identifier]).to eq(identifier)
        expect(Magma::Gnomon::Identifier.count).to eq(1)
        expect(Magma::Gnomon::Identifier.first[:identifier]).to eq(identifier)
      end

      context 'throws exception when' do
        it 'invalid rule name provided' do
          auth_header(:admin)
          post('/gnomon/labors/generate/alias/LABORS-LION-H2-C1')
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'non-admin tries to generate an identifier' do
          auth_header(:viewer)
          post('/gnomon/labors/generate/victim/LABORS-LION-H2-C1')
          expect(last_response.status).to eq(403)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'identifier already exists' do
          identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
          auth_header(:admin)
          post('/gnomon/labors/generate/victim/LABORS-LION-H2-C1')
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'identifier does not match rule' do
          auth_header(:admin)
          post('/gnomon/labors/generate/victim/LABORS-LION-H2-X1')
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end
      end
    end
  end
end
