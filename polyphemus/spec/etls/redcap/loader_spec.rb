describe Polyphemus::RedcapEtlScriptRunner do
  context 'entity iteration' do
    before(:each) do
      stub_magma_models(fixture: 'spec/fixtures/magma_test2_models.json')
      stub_redcap_test2_data
      @magma_client = Etna::Clients::Magma.new(host: MAGMA_HOST, token: TEST_TOKEN)
    end

    def run_loader(models)
      Redcap::Loader.new(
        {
          project_name: 'test2',
          models_to_build: models.keys.map(&:to_s),
          tokens: ["faketoken"],
          dateshift_salt: '123',
          redcap_host: REDCAP_HOST,
          magma_host: MAGMA_HOST,
          models: models
        }, 'test', @magma_client, STDERR
      ).run.first
    end

    it 'iterates across records' do
      Redcap::Model.define("Make").class_eval do
        def identifier(record_name)
          record_name
        end
      end

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

    it 'iterates across events' do
      Redcap::Model.define("Model").class_eval do
        def identifier(record_name, event_name)
          event_name
        end
      end

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

      Kernel.send(:remove_const,:Model)
    end

    it 'iterates across repeating instruments' do
      Redcap::Model.define("Year").class_eval do
        def identifier(record_name, event_name, (repeat_instrument, repeat_id) )
          [ record_name, event_name, repeat_id ].join(' ')
        end
      end

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

    it 'iterates across values' do
      Redcap::Model.define("Feature").class_eval do
        def patch(magma_record_name, record)
          record[:year] = magma_record_name.split('-').values_at(1,2,4).join(' ')
        end
      end

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
        def identifier(record_name, event_name, (repeat_instrument, repeat_id) )
          [ record_name, event_name, repeat_id ].join(' ')
        end
      end

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
