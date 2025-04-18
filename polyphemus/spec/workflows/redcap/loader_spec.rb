describe Redcap::Loader do
  def run_loader(models)
    Redcap::Loader.new(
      {
        project_name: 'test2',
        models_to_build: models.keys.map(&:to_s),
        tokens: ["faketoken"],
        redcap_host: REDCAP_HOST,
        magma_host: MAGMA_HOST,
        config: models
      }, 'test', @magma_client, STDERR
    ).run.first
  end

  context 'flat records' do
    before(:each) do
      stub_magma_models(fixture: 'spec/fixtures/magma_test2_models.json')
      @magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)

      Redcap::Model.define("Model").class_eval do
        def identifier(record_name, event_name, identifier_fields: nil)
          event_name
        end
      end
    end

    after(:each) do
      Kernel.send(:remove_const,:Model)
    end

    it 'uses flat record values to replace regular values' do
      stub_redcap(
        hash_including(content: 'metadata') => redcap_metadata(
          cars: [
            [ "record", "car_class" ]
          ]
        ),
        /eav/ => redcap_records(
          { field_name: "car_class" },
          [
            { "record": "Caudillac", "redcap_event_name": "El Corazon", "value": "1" },
            { "record": "Jatsun", "redcap_event_name": "Thunderer", "value": "2" },
          ]
        ).to_json,
        /flat/ => flat_records(
          redcap_records(
            { field_name: "car_class" },
            [
              { "record": "Caudillac", "redcap_event_name": "El Corazon", "value": "coupe" },
              { "record": "Jatsun", "redcap_event_name": "Thunderer", "value": "sedan" },
            ]
          )
        ).to_json
      )


      records = run_loader(
        model: {
          each: [ "record", "event" ],
          scripts: [
            {
              attributes: {
                type: "car_class"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:model])
      expect(records[:model].keys).to match_array(["Thunderer", "El Corazon"])
      expect(records[:model].values).to match_array([{type: "coupe"}, {type: "sedan"}])
    end

    it 'imposes an age cap for age configs' do
      stub_redcap(
        hash_including(content: 'metadata') => redcap_metadata(
          cars: [
            [ "record", "inventor_age" ]
          ]
        ),
        /eav/ => redcap_records(
          { field_name: "inventor_age" },
          [
            { "record": "Fjord", "redcap_event_name": "Model Eh", "value": "1" },
            { "record": "Verbie", "redcap_event_name": "Steamer", "value": "2" },
          ]
        ).to_json,
        /flat/ => flat_records(
          redcap_records(
            { field_name: "inventor_age" },
            [
              { "record": "Fjord", "redcap_event_name": "Model Eh", "value": "45" },
              { "record": "Verbie", "redcap_event_name": "Steamer", "value": "90" },
            ]
          )
        ).to_json
      )


      records = run_loader(
        model: {
          each: [ "record", "event" ],
          scripts: [
            {
              attributes: {
                inventor_age: {
                  "redcap_field": "inventor_age",
                  "value": "age"
                },
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:model])
      expect(records[:model].keys).to match_array(["Model Eh", "Steamer"])
      expect(records[:model].values).to match_array([{inventor_age: 45}, {inventor_age: 89}])
    end
  end


  context 'loading' do
    before(:each) do
      stub_magma_models(fixture: 'spec/fixtures/magma_test2_models.json')
      stub_redcap(
        hash_including(content: 'metadata') => redcap_metadata(
          cars: [
            [ "company_name" ],
            [ "dof", "date of founding" ],
            [ "car_class" ],
            [ "year", "calendar year", "date" ],
            {
              field_name: "feature",
              field_type: "checkbox",
              select_choices_or_calculations: redcap_choices("Window", "Door", "Spoiler", "Radio", "Engine", "Carbuerator", "Wheels")
            },
            [ "award_name", "Awards" ]
          ]
        )
      )
      @magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)
    end

    it 'combines records by record name' do
      Redcap::Model.define("Make").class_eval do
        def identifier(record_name, identifier_fields: nil)
          record_name
        end

        def offset_id(record_name)
          record_name
        end
      end

      stub_redcap(
        /fields/ => (redcap_records(
          { redcap_event_name: "Base", field_name: "dof" },
          [
            { record: "Jatsun", "value": "1956-02-03" },
            { "record": "ToyoT", "value": "1933-08-11" },
            { "record": "Caudillac", "value": "1901-03-01" },
          ]
        )+redcap_records(
          { redcap_event_name: "Base", field_name: "country" },
          [
            { record: "Jatsun", "value": "Japan" },
            { "record": "ToyoT", "value": "Japan" },
            { "record": "Caudillac", "value": "Spain" },
          ]
        )).to_json
      )

      records = run_loader(
        make: {
          each: [ "record" ],
          scripts: [
            {
              attributes: {
                date_of_founding: "dof"
              }
            },
            {
              attributes: {
                country: "country"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:make])
      expect(records[:make].keys).to match_array(["Jatsun", "ToyoT", "Caudillac"])
      expect(records[:make].values.map(&:keys)).to all( match_array([:country, :date_of_founding]) )

      Kernel.send(:remove_const,:Make)
    end
  end

  context 'entity iteration' do
    before(:each) do
      stub_magma_models(fixture: 'spec/fixtures/magma_test2_models.json')
      stub_redcap(
        hash_including(content: 'metadata') => redcap_metadata(
          cars: [
            [ "company_name" ],
            [ "dof", "date of founding" ],
            [ "car_class" ],
            [ "in_black" ],
            [ "year", "calendar year", "date" ],
            {
              field_name: "feature",
              field_type: "checkbox",
              select_choices_or_calculations: redcap_choices("Window", "Door", "Spoiler", "Radio", "Engine", "Carbuerator", "Wheels")
            },
            [ "award_name", "Awards" ]
          ]
        )
      )
      @magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)
    end

    it 'iterates across records' do
      Redcap::Model.define("Make").class_eval do
        def identifier(record_name, identifier_fields: nil)
          record_name
        end

        def offset_id(record_name)
          record_name
        end
      end

      stub_redcap(
        /fields/ => redcap_records(
          { redcap_event_name: "Base", field_name: "dof" },
          [
            { record: "Jatsun", "value": "1956-02-03" },
            { "record": "ToyoT", "value": "1933-08-11" },
            { "record": "Caudillac", "value": "1901-03-01" },
          ]
        ).to_json
      )

      records = run_loader(
        make: {
          each: [ "record" ],
          scripts: [
            {
              attributes: {
                date_of_founding: "dof"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:make])
      expect(records[:make].keys).to match_array(["Jatsun", "ToyoT", "Caudillac"])

      Kernel.send(:remove_const,:Make)
    end

    it 'correctly catch variants of Yes/No booleans' do
      Redcap::Model.define("Make").class_eval do
        def identifier(record_name, identifier_fields: nil)
          record_name
        end

        def offset_id(record_name)
          record_name
        end
      end

      stub_redcap(
        /fields/ => redcap_records(
          { redcap_event_name: "Base", field_name: "in_black" },
          [
            { record: "Jatsun", "value": "Yes" },
            { "record": "ToyoT", "value": "No" },
            { "record": "Caudillac", "value": "True" },
            { "record": "Hunda", "value": "False" },
          ]
        ).to_json
      )

      records = run_loader(
        make: {
          each: [ "record" ],
          scripts: [
            {
              attributes: {
                available_in_black: "in_black"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:make])
      expect(records[:make]["Jatsun"][:available_in_black]).to eq(true)
      expect(records[:make]["ToyoT"][:available_in_black]).to eq(false)
      expect(records[:make]["Caudillac"][:available_in_black]).to eq(true)
      expect(records[:make]["Hunda"][:available_in_black]).to eq(false)

      Kernel.send(:remove_const,:Make)
    end

    it 'iterates across events' do
      Redcap::Model.define("Model").class_eval do
        def identifier(record_name, event_name, identifier_fields: nil)
          event_name
        end
      end

      stub_redcap(
        /fields/ => redcap_records(
          { field_name: "car_class" },
          [
            { "record": "Caudillac", "redcap_event_name": "El Corazon", "value": "sedan" },
            { "record": "Jatsun", "redcap_event_name": "Thunderer", "value": "coupe" },
          ]
        ).to_json
      )

      records = run_loader(
        model: {
          each: [ "record", "event" ],
          scripts: [
            {
              attributes: {
                type: "car_class"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:model])
      expect(records[:model].keys).to match_array(["Thunderer", "El Corazon"])
      expect(records[:model].values).to match_array([{type: "coupe"}, {type: "sedan"}])

      Kernel.send(:remove_const,:Model)
    end

    it 'iterates across repeating instruments' do
      Redcap::Model.define("Year").class_eval do
        def identifier(record_name, event_name, (repeat_instrument, repeat_id), identifier_fields: nil )
          [ record_name, event_name, repeat_id ].join(' ')
        end
      end

      stub_redcap(
        /fields/ => redcap_records(
          { "redcap_repeat_instrument": "model_versions", "field_name": "year" },
          [
            {
              "record": "Jatsun", "redcap_event_name": "Thunderer",
              "redcap_repeat_instance": 1, "value": "1968"
            },
            {
              "record": "Jatsun", "redcap_event_name": "Thunderer",
              "redcap_repeat_instance": 2, "value": "1969"
            },
            {
              "record": "ToyoT", "redcap_event_name": "CorolloroC",
              "redcap_repeat_instance": 1, "value": "1979"
            },
            {
              "record": "ToyoT", "redcap_event_name": "CorolloroC",
              "redcap_repeat_instance": 2, "value": "1980"
            }
          ]
        ).to_json
      )

      records = run_loader(
        year: {
          each: [ "record", "event", "repeat" ],
          scripts: [
            {
              attributes: {
                calendar_year: "year"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:year])
      expect(records[:year].keys).to match_array(["Jatsun Thunderer 1", "Jatsun Thunderer 2", "ToyoT CorolloroC 1", "ToyoT CorolloroC 2"])

      Kernel.send(:remove_const,:Year)
    end

    it 'correctly handles repeat_instance id for eav and flat across repeating instruments with identifier_fields' do
      Redcap::Model.define("Award").class_eval do
        def identifier(record_name, event_name = nil, identifier_fields: nil )
          [
            "::temp", identifier_fields[:special_id], rand(36 ** 8).to_s(36),
          ].compact.join("-")
        end

        def patch(id, record)
          record[:model] = id.split('-')[1]
        end
      end

      stub_redcap(
        /eav/ => redcap_records(
          { "redcap_repeat_instrument": "awards" },
          [
            {
              "record": "Jatsun", "field_name": "name",
              "redcap_repeat_instance": 1, "value": "Gold"
            },
            {
              "record": "Jatsun", "field_name": "name",
              "redcap_repeat_instance": '2', "value": "Silver"
            },
            {
              "record": "ToyoT", "field_name": "name",
              "redcap_repeat_instance": 1, "value": "Bronze"
            },
            {
              "record": "ToyoT", "field_name": "name",
              "redcap_repeat_instance": '2', "value": "Platinum"
            }
          ]
        ).to_json,
        /flat/ => redcap_records(
          { },
          [
            {
              "redcap_repeat_instrument": "",
              "record_id": "Jatsun", "redcap_repeat_instance": "",
              "special_id": "brand1", "company_name": "Jatsun"
            },
            {
              "redcap_repeat_instrument": "awards",
              "record_id": "Jatsun", "name": "Gold",
              "redcap_repeat_instance": 1, "company_name": "Jatsun",
              "special_id": ""
            },
            {
              "redcap_repeat_instrument": "awards",
              "record_id": "Jatsun", "name": "Silver",
              "redcap_repeat_instance": 2, "company_name": "Jatsun",
              "special_id": ""
            },
            {
              "redcap_repeat_instrument": "awards",
              "record_id": "ToyoT", "name": "Bronze",
              "redcap_repeat_instance": 1, "company_name": "ToyoT",
              "special_id": ""
            },
            {
              "redcap_repeat_instrument": "awards",
              "record_id": "ToyoT", "name": "Platinum",
              "redcap_repeat_instance": 2, "company_name": "ToyoT",
              "special_id": ""
            },
            {
              "redcap_repeat_instrument": "",
              "record_id": "ToyoT", "redcap_repeat_instance": "",
              "special_id": "brand2", "company_name": "ToyoT"
            },
          ]
        ).to_json
      )

      records = run_loader(
        award: {
          each: [ "record", "repeat" ],
          identifier_fields: ["special_id"],
          scripts: [
            {
              attributes: {
                name: "name"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:award, :model])
      expect(records[:award].keys.length).to eq(4)
      expect(records[:model]["brand1"]["award"].length).to eq(2)
      expect(records[:model]["brand2"]["award"].length).to eq(2)

      Kernel.send(:remove_const,:Award)
    end

    it 'iterates across values' do
      Redcap::Model.define("Feature").class_eval do
        def patch(magma_record_name, record)
          record[:year] = magma_record_name.split('-').values_at(1,2,4).join(' ')
        end
      end

      stub_redcap(
        /fields/ => (redcap_records(
          { "redcap_repeat_instrument": "model_versions", "field_name": "feature", "redcap_event_name": "CorolloroC" },
          [
            { "record": "ToyoT", "redcap_repeat_instance": 1, "value": "Door" },
            { "record": "ToyoT", "redcap_repeat_instance": 1, "value": "Engine" },
            { "record": "ToyoT", "redcap_repeat_instance": 2, "value": "Wheels" }
          ]
        ) + redcap_records(
          { "redcap_repeat_instrument": "model_versions", "field_name": "year", "record": "ToyoT", "redcap_event_name": "CorolloroC" },
          [
            { "redcap_repeat_instance": 1, "value": "1979" },
            { "redcap_repeat_instance": 2, "value": "1980" }
          ]
        )).to_json
      )

      records = run_loader(
        feature: {
          each: [ "record", "event", "repeat", { "field" => "feature" } ],
          scripts: [
            {
              attributes: {
                name: {
                  "redcap_field": "feature",
                  "value": "value"
                },
                value: "year"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:feature, :year])
      expect(records[:feature].values).to match_array([
        {:name=>"Door", :value=>"1979", :year=>"ToyoT CorolloroC 1"},
        {:name=>"Engine", :value=>"1979", :year=>"ToyoT CorolloroC 1"},
        {:name=>"Wheels", :value=>"1980", :year=>"ToyoT CorolloroC 2"}
      ])

      Kernel.send(:remove_const, :Feature)
    end

    it 'combines across values' do
      Redcap::Model.define("Year").class_eval do
        def identifier(record_name, event_name, (repeat_instrument, repeat_id), identifier_fields: nil )
          [ record_name, event_name, repeat_id ].join(' ')
        end
      end

      stub_redcap(
        /fields/ => (redcap_records(
          { "redcap_repeat_instrument": "model_versions", "field_name": "feature", "redcap_event_name": "CorolloroC" },
          [
            { "record": "ToyoT", "redcap_repeat_instance": 1, "value": "Door" },
            { "record": "ToyoT", "redcap_repeat_instance": 1, "value": "Engine" },
            { "record": "ToyoT", "redcap_repeat_instance": 2, "value": "Wheels" }
          ]
        ) + redcap_records(
          { "redcap_repeat_instrument": "model_versions", "field_name": "year", "record": "ToyoT", "redcap_event_name": "CorolloroC" },
          [
            { "redcap_repeat_instance": 1, "value": "1979" },
            { "redcap_repeat_instance": 2, "value": "1980" }
          ]
        ) + redcap_records(
          { "record": "Jatsun", "redcap_event_name": "Thunderer", "redcap_repeat_instrument": "model_versions", "field_name": "year" },
          [
            { "redcap_repeat_instance": 1, "value": "1968" },
            { "redcap_repeat_instance": 2, "value": "1969" }
          ]
        ) ).to_json
      )

      records = run_loader(
        year: {
          each: [ "record", "event", "repeat" ],
          scripts: [
            {
              attributes: {
                calendar_year: "year",
                feature_list: {
                  "redcap_field": "feature",
                  "value": "combine",
                  "combine": ", "
                }
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:year])
      expect(records[:year].values).to match_array(
        [
          {:calendar_year=>"1968"},
          {:calendar_year=>"1969"},
          {:calendar_year=>"1979", :feature_list=>"Door, Engine"},
          {:calendar_year=>"1980", :feature_list=>"Wheels"}
        ]
      )

      Kernel.send(:remove_const, :Year)
    end

    it 'iterates with restrictions' do
      Redcap::Model.define("Award").class_eval do
        def patch(magma_record_name, record)
          record[:model] = magma_record_name.split('-')[2]
        end
      end

      stub_redcap(
        /fields/ => (redcap_records(
          { "record": "ToyoT", "redcap_event_name": "CorolloroC", "redcap_repeat_instance": 1, "field_name": "award_name" },
          [
            { "redcap_repeat_instrument": "awards", "value": "Car of the Year 1977" },
            { "redcap_repeat_instrument": "rejected", "value": "Car of the Year 1978" }
          ]
        ) ).to_json
      )

      records = run_loader(
        award: {
          each: [ "record", "event", { repeat: /awards/  } ],
          scripts: [
            {
              attributes: {
                name: "award_name"
              }
            }
          ]
        }
      )

      expect(records.keys).to match_array([:award, :model])
      expect(records[:award].values).to match_array([
        {model: "CorolloroC", name: "Car of the Year 1977"}
      ])
      Kernel.send(:remove_const,:Award)
    end
  end
end
