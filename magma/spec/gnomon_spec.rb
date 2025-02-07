
require 'json'


def create_grammar(params={})
  grammar = create(:grammar, { project_name: 'labors', version_number: 1, config: {}, comment: 'update' }.merge(params))
end

describe GnomonController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it 'complains if there is no grammar' do
    auth_header(:viewer)
    get('/gnomon/labors')

    expect(last_response.status).to eq(422)
  end

  it 'gets the most recent grammar' do
    grammar = create_grammar(version_number: 1, config: {})
    grammar2 = create_grammar(version_number: 2, config: VALID_GRAMMAR_CONFIG)
    auth_header(:viewer)
    get('/gnomon/labors')

    expect(last_response.status).to eq(200)
    expect(json_body.to_json).to eq(grammar2.to_hash.to_json)
  end

  it 'sets a new grammar' do
    grammar = create_grammar

    config = VALID_GRAMMAR_CONFIG
    auth_header(:admin)
    json_post('/gnomon/labors', config: config, comment: 'eh')

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
    json_post('/gnomon/labors', config: config, comment: 'eh')

    expect(last_response.status).to eq(422)
    expect(json_body).to match_array(errors: [
      "root is missing required keys: tokens, rules",
      "property '/text' is invalid: error_type=schema",
      "No separator token defined!"
    ])
    expect(Magma::Gnomon::Grammar.count).to eq(0)
  end

  it 'decomposes an identifier' do
    Timecop.freeze
    grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
    identifier = create_identifier("The Twelve Labors of Hercules", rule: 'project', grammar: grammar)
    identifier2 = create_identifier("The Nemean Lion", rule: 'labor', grammar: grammar)
    record = create(:project, name: "The Twelve Labors of Hercules")
    auth_header(:viewer)
    get('/gnomon/labors/decompose/LABORS-LION-H2-C1')

    expect(last_response.status).to eq(200)
    expect(json_body).to eq(
      rule_name: "victim",
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
    grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
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
        @grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
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

  context 'rule API' do
    it 'throws exception when no grammar for project' do
      auth_header(:viewer)
      get('/gnomon/labors/rule/victim')
      expect(last_response.status).to eq(422)
    end

    context 'with grammar' do
      before(:each) do
        @grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
      end

      it 'returns a tokenization of the rule' do
        identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
        auth_header(:viewer)
        get('/gnomon/labors/rule/village')
        expect(last_response.status).to eq(200)
        expect(json_body[:rule]).to eq(
          [
            {label: "project", values: {LABORS: "The Twelve Labors of Hercules"}, name: "PROJ"},
            {label: "Separator", values: {'-': "# Separator"}, name: "SEP"},
            {label: "labor", values: {LION: "The Nemean Lion", HYDRA: "The Lernean Hydra"}, name: "LAB"},
            {label: "Separator", values: {'-': "# Separator"}, name: "SEP"},
            {label: "Village type", values: {V: "Village", H: "Hamlet"}, name: "VILL"},
            {name: "n", label: "village_counter"}
          ]
        )
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
        @grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
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
        @grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
      end

      it 'creates the identifier' do
        identifier = "LABORS-LION-H2-C1"
        expect(Magma::Gnomon::Identifier.count).to eq(0)
        auth_header(:admin)
        post("/gnomon/labors/generate/victim/#{identifier}")

        rules = {
          labor: "The Nemean Lion",
          project: "The Twelve Labors of Hercules",
          victim: "LABORS-LION-H2-C1",
          village: "LABORS-LION-H2"
        }
        expect(last_response.status).to eq(200)
        expect(json_body[:tokens].map(&:last).join).to eq("LABORS-LION-H2-C1")
        expect(json_body[:rules].to_h{|k,v| [ k,v[:name] ]}).to eq(rules)
        expect(Magma::Gnomon::Identifier.count).to eq(4)
        expect(Magma::Gnomon::Identifier.select_map(
          [ :rule, :identifier ]
        ).to_h.symbolize_keys).to eq(rules)
      end

      it 'ignores identifiers that already exists' do
        identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
        auth_header(:admin)
        post('/gnomon/labors/generate/victim/LABORS-LION-H2-C1')
        expect(last_response.status).to eq(200)
        expect(json_body[:rules].size).to eq(4)
        expect(Magma::Gnomon::Identifier.count).to eq(4)
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

        it 'identifier does not match rule' do
          auth_header(:admin)
          post('/gnomon/labors/generate/victim/LABORS-LION-H2-X1')
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end
      end
    end
  end

  context 'bulk generate API' do

    it 'throws exception when no grammar for project' do
      auth_header(:admin)
      json_post('/gnomon/labors/generate', names: [])
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('No grammar found for project labors.')
    end

    context 'with grammar' do
      before(:each) do
        @grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)
      end

      it 'creates identifiers' do
        auth_header(:admin)
        json_post("/gnomon/labors/generate", names: [
          {
            rule_name: 'village',
            name: 'LABORS-HYDRA-H2'
          },
          {
            rule_name: 'victim',
            name: 'LABORS-LION-H2-C1'
          }
        ])

        expect(last_response.status).to eq(200)
        expect(Magma::Gnomon::Identifier.count).to eq(6)
        expect(Magma::Gnomon::Identifier.select_map(
          [ :rule, :identifier ]
        )).to match_array([
          ["labor", "The Lernean Hydra"],
          ["labor", "The Nemean Lion"],
          ["project", "The Twelve Labors of Hercules"],
          ["victim", "LABORS-LION-H2-C1"],
          ["village", "LABORS-HYDRA-H2"],
          ["village", "LABORS-LION-H2"]
        ])
      end

      it 'ignores identifiers that already exists' do
        identifier = create_identifier("LABORS-LION-H2-C1", rule: 'victim', grammar: @grammar)
        auth_header(:admin)
        names = [{rule_name: 'victim', name: 'LABORS-LION-H2-C1'}]
        json_post('/gnomon/labors/generate', names: names)
        expect(last_response.status).to eq(200)
        expect(json_body[:existing]).to eq(names)
        expect(Magma::Gnomon::Identifier.count).to eq(4)
      end


      context 'throws exception when' do
        it 'invalid rule name provided' do
          auth_header(:admin)
          json_post('/gnomon/labors/generate', names: [{ rule_name: 'alias', name: 'LABORS-LION-H2-C1'}])
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'non-admin tries to generate an identifier' do
          auth_header(:viewer)
          json_post('/gnomon/labors/generate', names: [{ rule_name: 'victim', name: 'LABORS-LION-H2-C1'}])
          expect(last_response.status).to eq(403)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end

        it 'identifier does not match rule' do
          auth_header(:admin)
          json_post('/gnomon/labors/generate', names: [{ rule_name: 'victim', name: 'LABORS-LION-H2-X1'}])
          expect(last_response.status).to eq(422)
          expect(Magma::Gnomon::Identifier.count).to eq(0)
        end
      end
    end
  end

  context 'rules API' do
    before(:each) do
      grammar = create_grammar(version_number: 1, config: {})
      grammar2 = create_grammar(version_number: 2, config: VALID_GRAMMAR_CONFIG)
      grammar3 = create_grammar(project_name: 'toils', version_number: 1, config: VALID_GRAMMAR_CONFIG)
    end

    it 'lists all of the rules for each requested project' do
      auth_header(:superuser)
      post('/gnomon/rules', project_names: [ 'labors', 'toils' ])
      expect(last_response.status).to eq(200)

      expect(json_body[:rules]).to eq(
       labors: {
         labor: "^(The Nemean Lion|The Lernean Hydra)$",
         project: "^(The Twelve Labors of Hercules)$",
         victim: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)(-)(S|C)(\\d+)$",
         village: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)$"
       },
       toils: {
         labor: "^(The Nemean Lion|The Lernean Hydra)$",
         project: "^(The Twelve Labors of Hercules)$",
         victim: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)(-)(S|C)(\\d+)$",
         village: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)$"
       }
      )
    end

    it 'lists all of the rules for a given project' do
      auth_header(:viewer)
      get('/gnomon/labors/rules')
      expect(last_response.status).to eq(200)

      expect(json_body[:rules]).to eq(
       labor: "^(The Nemean Lion|The Lernean Hydra)$",
       project: "^(The Twelve Labors of Hercules)$",
       victim: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)(-)(S|C)(\\d+)$",
       village: "^(LABORS)(-)(LION|HYDRA)(-)(V|H)(\\d+)$"
      )
    end
  end
  end
end

describe Magma::Gnomon::Identifier do
  context 'backfills' do
    it 'only good identifiers' do
      grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)

      victim1 = create(:victim, name: 'Outis Koutsonadis')
      victim2 = create(:victim, name: 'Susan Doe')
      victim3 = create(:victim, name: 'LABORS-LION-H2-C1')
      victim4 = create(:victim, name: 'LABORS-HYDRA-H1-C1')

      expect {
        Magma::Gnomon::Identifier.backfill(
          'labors', 'victim', 'eurystheus@twelve-labors.org|Eurystheus',
          dry_run: false
        )

      }.to output("2 ids backfilled, 2 do not match current grammar.\nBad Identifiers:\nOutis Koutsonadis, Susan Doe\n").to_stdout

      expect(Magma::Gnomon::Identifier.select_map(:identifier)).to match_array([
        'LABORS-LION-H2-C1', 'LABORS-HYDRA-H1-C1'
      ])
    end

    it 'dry-runs' do
      grammar = create_grammar(config: VALID_GRAMMAR_CONFIG)

      victim1 = create(:victim, name: 'Outis Koutsonadis')
      victim2 = create(:victim, name: 'Susan Doe')
      victim3 = create(:victim, name: 'LABORS-LION-H2-C1')
      victim4 = create(:victim, name: 'LABORS-HYDRA-H1-C1')

      expect {
        Magma::Gnomon::Identifier.backfill(
          'labors', 'victim', 'eurystheus@twelve-labors.org|Eurystheus'
        )

      }.to output("2 ids backfilled, 2 do not match current grammar.\nBad Identifiers:\nOutis Koutsonadis, Susan Doe\n").to_stdout

      expect(Magma::Gnomon::Identifier.count).to eq(0)
    end
  end
end
