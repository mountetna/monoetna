describe Polyphemus::RedcapEtlScriptRunner do
  TEST_REDCAP_CONFIG = {
    model_one: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            birthday: "date_of_birth",
            graduation_date: "commencement",
            name: "name"
          }
        }
      ]
    },
    model_two: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            label: {
              redcap_field: "today",
              value: "label"
            },
            yesterday: {
              redcap_field: "today",
              value: "value"
            }
          }
        }
      ]
    },
    stats: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            height: {
              redcap_field: "height",
              value: "value"
            },
            weight: {
              redcap_field: "weight",
              value: "value"
            }
          }
        }
      ]
    },
    citation: {
      each: [ "record", "event" ],
      invert: true,
      scripts: [
        {
          attributes: {
            name: {
              match: 'citation-(111|abc)',
              value: "none"
            },
            date: {
              redcap_field: "citation_date",
              value: "value"
            }
          }
        }
      ]
    }
  }

  NO_OFFSET_ID_REDCAP_CONFIG = {
    bad_model: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            birthday: "date_of_birth",
          }
        }
      ]
    }
  }

  ALTERNATE_ID_REDCAP_CONFIG = {
    model_with_alternate_id: {
      each: [ "record" ],
      identifier_fields: [ "date_of_birth" ],
      scripts: [
        {
          attributes: {
            name: "name",
          }
        }
      ]
    }
  }

  CONFIG_WITH_FILTERS = {
    model_one: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            birthday: "date_of_birth",
            graduation_date: "commencement",
            name: "name"
          },
          filters: [{
            redcap_field: "rain_date",
            exists: true
          }]
        }
      ]
    },
    model_two: {
      each: [ "record" ],
      scripts: [
        {
          attributes: {
            label: {
              redcap_field: "today",
              value: "label"
            },
            yesterday: {
              redcap_field: "today",
              value: "value"
            }
          },
          filters: [{
            redcap_field: "some_day",
            equals: "2022-01-01"
          }]
        }
      ]
    },
  }

  context 'dateshifts' do
    before do
      stub_magma_models
      copy_redcap_project
    end

    it 'throws exception if salt not provided' do
      expect {
        Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "all",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: nil,
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          config: TEST_REDCAP_CONFIG
        )
      }.to raise_error(RuntimeError, "No dateshift_salt provided, please provide one.")
    end

    it 'uses the provided salt' do
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "all",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: TEST_REDCAP_CONFIG
      )

      system_config = redcap_etl.system_config

      expect(system_config[:dateshift_salt]).to eq('123')
    end

    it 'throws exception if model with date_time attribute does not define offset_id' do
      stub_redcap_data(:essential_data)
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "all",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: NO_OFFSET_ID_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      expect {
        redcap_etl.run(magma_client: magma_client)
      }.to raise_error(RuntimeError, "offset_id() needs to be implemented for the test project, bad_model class. It should return the patient / subject identifier.")
    end
  end

  context 'redcap tokens' do
    before do
      stub_magma_models
      stub_redcap_data
      copy_redcap_project
    end

    it 'throws exception if not provided' do
      expect {
        Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "all",
          redcap_tokens: [],
          dateshift_salt: "123",
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          config: TEST_REDCAP_CONFIG
        )
      }.to raise_error(RuntimeError, "Must provide at least one REDCap token.")
    end
  end

  context 'loads' do
    before do
      stub_magma_update_dry_run
      stub_magma_models
      stub_redcap_data
      copy_redcap_project
    end

    let(:raw_data) {
      arry1 = JSON.parse(File.read('spec/fixtures/redcap_mock_data_calendar.json'), symbolize_names: true)
      arry2 = JSON.parse(File.read('spec/fixtures/redcap_mock_data_essential_data.json'), symbolize_names: true)

      return arry1 + arry2
    }

    let(:raw_template) {
      JSON.parse(File.read('spec/fixtures/redcap_mock_metadata.json'), symbolize_names: true)
    }

    def data_for_id(id)
      raw_data.each do |record|
        return record if record[:record] == id
      end
      return nil
    end

    def template_for_field(field)
      raw_template.each do |record|
        return record if record[:field_name] == field
      end
      return nil
    end

    it 'all models' do
      stub_redcap_data(:essential_data)
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "all",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: TEST_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      raw_data_456 = data_for_id('456')
      raw_data_000 = data_for_id('000')
      raw_data_123 = data_for_id('123')
      raw_data_321 = data_for_id('321')
      id_456 = temp_id(records, "456")
      id_000 = temp_id(records, "000")
      id_123 = "123"
      id_321 = "321"

      # Dateshift checks
      expect(records[:model_one][id_456][:birthday]).not_to eq(raw_data_456[:value])
      expect(records[:model_one][id_456][:birthday].start_with?('1899')).to eq(true)
      expect(records[:model_one][id_000][:graduation_date]).not_to eq(raw_data_000[:value])
      expect(records[:model_one][id_000][:graduation_date].start_with?('2021')).to eq(true)
      expect(records[:model_two][id_123][:yesterday]).not_to eq(raw_data_123[:value])
      expect(records[:model_two][id_123][:yesterday].start_with?('2019')).to eq(true)
      expect(records[:model_two][id_321][:yesterday]).not_to eq(raw_data_321[:value])
      expect(records[:model_two][id_321][:yesterday].start_with?('2020')).to eq(true)

      # ensure containing records works for model_two
      injected_parent_id = "#{id_123}-one"
      expect(records[:model_one][injected_parent_id][:parent_model]).to eq("#{injected_parent_id}-parent")
    end

    it 'specific models' do
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "model_one",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: TEST_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      expect(records.keys.include?(:model_one)).to eq(true)
      expect(records[:model_one].keys.length).to eq(6)
      expect(records.keys.include?(:model_two)).to eq(false)
    end

    it 'form label attributes' do
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "model_two",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: TEST_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      raw_template_today = template_for_field('today')
      id_123 = "123"

      expect(records[:model_two][id_123][:label]).to eq(raw_template_today[:field_label])
    end

    it 'multiple REDCap projects' do
      stub_redcap_multi_project_records

      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "model_one",
        redcap_tokens: "#{REDCAP_TOKEN},#{REDCAP_TOKEN.reverse}",
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: TEST_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      expect(records.keys.include?(:model_one)).to eq(true)
      expect(records.keys.include?(:model_two)).to eq(false)

      expect(records[:model_one].keys.length).to eq(6)
    end

    it 'when using alternate REDCap fields for ids' do
      stub_redcap_data(:essential_data)
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "model_with_alternate_id",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: ALTERNATE_ID_REDCAP_CONFIG
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      expect(records.keys.include?(:model_with_alternate_id)).to eq(true)

      # Only finds records with all the fields in the CONFIG, which
      #   should be :date_of_birth
      expect(records[:model_with_alternate_id].keys.length).to eq(2)

      data_789 = data_for_id('789')
      data_987 = data_for_id('987')

      expected_id_789 = temp_id(records, data_789[:date_of_birth])
      expected_id_987 = temp_id(records, data_987[:date_of_birth])

      expect(records[:model_with_alternate_id].key?(expected_id_789)).to eq(true)
      expect(records[:model_with_alternate_id].key?(expected_id_987)).to eq(true)
    end

    it "filters records based on non-mapped redcap attribute" do
      redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
        project_name: 'test',
        model_names: "all",
        redcap_tokens: REDCAP_TOKEN,
        dateshift_salt: '123',
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: CONFIG_WITH_FILTERS
      )

      magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      records = redcap_etl.run(magma_client: magma_client)

      raw_data_000 = data_for_id('000')
      raw_data_321 = data_for_id('321')
      id_456 = temp_id(records, "456")
      id_000 = temp_id(records, "000")
      id_123 = "123"
      id_321 = "321"

      # Verify that only record 000 for model one passes the filter
      expect(records[:model_one].keys.include?(id_456)).to eq(false)
      expect(records[:model_one][id_000][:graduation_date]).not_to eq(raw_data_000[:value])
      expect(records[:model_one][id_000][:graduation_date].start_with?('2021')).to eq(true)
      expect(records[:model_one][id_000].keys.include?(:rain_date)).to eq(false) # filter attribute

      # Verify that only record 321 for model two passes the filter, and that
      # ensure containing records works for model_two still
      expect(records[:model_two].keys.include?(id_123)).to eq(false)
      expect(records[:model_two][id_321][:yesterday]).not_to eq(raw_data_321[:value])
      expect(records[:model_two][id_321][:yesterday].start_with?('2020')).to eq(true)
      expect(records[:model_two][id_321].keys.include?(:some_day)).to eq(false) # filter attribute
      injected_parent_id = "#{id_321}-one"
      expect(records[:model_one][injected_parent_id][:parent_model]).to eq("#{injected_parent_id}-parent")
    
    end

    context("mode == nil") do
      it 'updates records even if not in Magma with regular model' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "all",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        expect(records.keys.include?(:model_one)).to eq(true)
        expect(records.keys.include?(:model_two)).to eq(true)

        expect(records[:model_two].keys).to match_array(["123", "321", "abc"])

        # Null value keys do not appear
        id_456 = temp_id(records, "456")
        expect(records[:model_one][id_456].keys.include?(:graduation_date)).to eq(false)
      end

      it 'updates records even if not in Magma with table model' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "stats",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)
        id_abc = temp_id(records, "abc")

        expect(records.keys.include?(:stats)).to eq(true)
        expect(records[:stats].keys.length).to eq(1)
        expect(records[:stats].keys).to match_array([id_abc])

        expect(records.keys.include?(:model_two)).to eq(true)
        expect(records[:model_two].keys).to match_array(["abc"])
        expect(records[:model_two]["abc"]["stats"]).to eq(records[:stats].keys)
      end
    end

    context("mode == existing") do
      it 'updates only records in Magma with regular model' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "model_two",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: "existing",
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        # Even though the parent record is attempted to be made,
        #   it should be removed because it's not an "existing" magma
        #   record per the fixture.
        expect(records[:model_one].keys.length).to eq(0)
        expect(records.keys.include?(:model_two)).to eq(true)
        expect(records[:model_two].keys).to eq(["123"])
      end

      it 'updates only records in Magma with table model' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "stats",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: 'existing',
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        expect(records.keys.include?(:stats)).to eq(true)
        expect(records[:stats].keys.length).to eq(0)
      end


    end

    context("mode == strict") do
      it 'shows null values on regular models' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "model_one",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: "strict",
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        expect(records.keys.include?(:model_one)).to eq(true)
        # Make sure null values appear now
        id_456 = temp_id(records, "456")
        expect(records[:model_one][id_456].keys.include?(:graduation_date)).to eq(true)
        expect(records[:model_one][id_456][:graduation_date]).to eq(nil)
      end

      it 'blanks parent records in Magma that are not in REDCap when updating table model' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "stats",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: "strict",
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        expect(records.keys.include?(:model_two)).to eq(true)
        expect(records[:model_two].keys.length).to eq(2)
        expect(records[:model_two]["abc"]["stats"]).not_to eq([])
        expect(records[:model_two]["123"]["stats"]).to eq([])
      end

      it 'blanks "removed" magma records in strict mode' do
        # Not a great test ... can't figure out how to test or mock for
        #   a process spun out in a different Thread.
        stub_magma_update_json

        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "stats",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: "strict",
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        id_abc = temp_id(records, "abc")

        expect(records.keys).to match_array([:stats, :model_two])
        expect(records[:model_two].keys).to match_array(['123', 'abc'])
        expect(records[:stats].keys).to match_array([id_abc])

        expect(records[:model_two]['abc']['stats']).not_to eq([])
        expect(records[:model_two]['123']['stats']).to eq([])
      end
    end

    context 'invert loading' do
      it 'invert loads redcap data into existing magma records' do
        redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
          project_name: 'test',
          model_names: "citation",
          redcap_tokens: REDCAP_TOKEN,
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          mode: "strict",
          config: TEST_REDCAP_CONFIG
        )

        magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

        records = redcap_etl.run(magma_client: magma_client)

        expect(records.keys).to match_array([:citation])
        expect(records[:citation].keys).to match_array(["citation-abc-1", "citation-abc-2", "citation-111-1"])
        expect(records[:citation].values.map{|r| r[:date]}).to all( match(/20..-..-../))
      end
    end
  end
end
