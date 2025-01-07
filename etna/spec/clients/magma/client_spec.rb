describe Etna::Clients::Magma do
  around(:each) do |example|
    VCR.turned_off { example.run }
  end

  before(:each) do
    WebMock.disable_net_connect!
  end

  let(:update) do
    {
      project_name: 'labors',
      revisions: {
        monster: {
          'Nemean Lion': {
            species: 'lion'
          },
          'Lernean Hydra': {
            species: 'hydra'
          },
          'Ceryneian Hind': {
            species: 'deer'
          }
        }
      }
    }
  end

  let(:magma_client) do
    Etna::Clients::Magma.new(
        host: 'https://magma.test',
        token: TEST_TOKEN
    )
  end

  let(:monster_model) do
    {
      models: {
        "monster": {
          "documents": {},
          "template": {
            "parent": "labor",
            "attributes": {
              "species": {
                "name": "species",
                "attribute_name": "species",
                "restricted": false,
                "read_only": false,
                "hidden": false,
                "validation": nil,
                "attribute_type": "string"
              },
              "victim": {
                "name": "victim",
                "attribute_name": "victim",
                "restricted": false,
                "read_only": false,
                "hidden": false,
                "validation": nil,
                "attribute_type": "table",
                "link_model_name": "victim",
                "link_attribute_name": "monster"
              }
            }
          }
        },
        "victim": {
          "documents": {},
          "template": {
            "parent": "monster",
            "attributes": {
              "monster": {
                "name": "monster",
                "attribute_name": "monster",
                "restricted": false,
                "read_only": false,
                "hidden": false,
                "validation": nil,
                "attribute_type": "parent",
                "link_model_name": "monster",
                "link_attribute_name": "victim"
              },
              "age": {
                "name": "age",
                "attribute_name": "age",
                "restricted": false,
                "read_only": false,
                "hidden": false,
                "validation": nil,
                "attribute_type": "integer"
              },
              "village": {
                "name": "village",
                "attribute_name": "village",
                "restricted": false,
                "read_only": false,
                "hidden": false,
                "validation": nil,
                "attribute_type": "string"
              }
            }
          }
        }
      }
    }
  end

  it 'posts an update' do
    stub_request(:post, /#{MAGMA_HOST}\/update$/).to_return(body: '{}')
    response = magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(update))

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including(update))
  end

  context 'with paged updates' do
    context 'for non-table models' do
      let(:response1) do
        {
          models: {
            monster: {
              documents: update[:revisions][:monster].slice(:'Nemean Lion', :'Lernean Hydra')
            }
          }
        }
      end

      let(:response2) do
        {
          models: {
            monster: {
              documents: update[:revisions][:monster].slice(:'Ceryneian Hind')
            }
          }
        }
      end

      it 'posts a paged update' do
        stub_request(:post, /#{MAGMA_HOST}\/update$/).to_return([
          { body: response1.to_json }, {body: response2.to_json}
        ])
        stub_request(:post, /#{MAGMA_HOST}\/retrieve$/).to_return(body: monster_model.to_json)
        response = magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(update), 2)

        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).twice
        expect(response.models.model("monster").documents.document_keys).to eq(["Nemean Lion", "Lernean Hydra", "Ceryneian Hind"])
      end
    end

    context 'for table models' do
      let(:blank_table_update) do
        {
          project_name: 'labors',
          revisions: {
            'monster' => {
              'Nemean Lion' => {
                'victim' => [ '::temp-victim1', '::temp-victim2', '::temp-victim3', '::temp-victim4', '::temp-victim5' ]
              },
              'Lernean Hydra' => {
                'victim' => [ '::temp-victim6', '::temp-victim7', '::temp-victim8', '::temp-victim9' ]
              }
            },
            'victim' => {
              '::temp-victim1' => { monster: 'Nemean Lion', age: '92', village: 'Nemea' },
              '::temp-victim2' => { monster: 'Nemean Lion', age: '26', village: 'Cleonae' },
              '::temp-victim3' => { monster: 'Nemean Lion', age: '25', village: 'Cleonae' },
              '::temp-victim4' => { monster: 'Nemean Lion', age: '81', village: 'Nemea' },
              '::temp-victim5' => { monster: 'Nemean Lion', age: '1', village: 'Nemea' },
              '::temp-victim6' => { monster: 'Lernean Hydra', age: '20', village: 'Lerna' },
              '::temp-victim7' => { monster: 'Lernean Hydra', age: '22', village: 'Lerna' },
              '::temp-victim8' => { monster: 'Lernean Hydra', age: '42', village: 'Lerna' },
              '::temp-victim9' => { monster: 'Lernean Hydra', age: '60', village: 'Lerna' }
            }
          }
        }
      end

      it 'posts a paged update appending to tables' do
        revisions1 = {
          "monster" => blank_table_update[:revisions]["monster"].slice('Nemean Lion'),
          "victim" => blank_table_update[:revisions]["victim"].slice(
            "::temp-victim1", "::temp-victim2", "::temp-victim3", "::temp-victim4", "::temp-victim5"
          )
        }
        revisions2 = {
          "monster" => blank_table_update[:revisions]["monster"].slice('Lernean Hydra'),
          "victim" => blank_table_update[:revisions]["victim"].slice(
            "::temp-victim6", "::temp-victim7", "::temp-victim8", "::temp-victim9"
          )
        }
        stub_request(:post, /#{MAGMA_HOST}\/update$/).with(
          body: hash_including("revisions" => revisions1)
        ).to_return(body: {
          models: {
            monster: { documents: revisions1["monster"] },
            victim: { documents: revisions1["victim"] }
          }
        }.to_json)

        stub_request(:post, /#{MAGMA_HOST}\/update$/).with(
          body: hash_including("revisions" => revisions2)
        ).to_return(body: {
          models: {
            monster: { documents: revisions2["monster"] },
            victim: { documents: revisions2["victim"] }
          }
        }.to_json)

        stub_request(:post, /#{MAGMA_HOST}\/retrieve$/).to_return(body: monster_model.to_json)
        response = magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(blank_table_update), 2)

        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).twice
        expect(response.models.model("monster").documents.document_keys).to eq(["Nemean Lion", "Lernean Hydra", "Ceryneian Hind"])
      end

      it 'posts a paged update blanking tables' do
        stub_request(:post, /#{MAGMA_HOST}\/update$/).to_return([
          { body: response1.to_json }, {body: response2.to_json}
        ])
        stub_request(:post, /#{MAGMA_HOST}\/retrieve$/).to_return(body: monster_model.to_json)
        response = magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(update), 2)

        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).twice
        expect(response.models.model("monster").documents.document_keys).to eq(["Nemean Lion", "Lernean Hydra", "Ceryneian Hind"])
      end
    end
  end
end
