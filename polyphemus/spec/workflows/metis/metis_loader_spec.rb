describe Metis::Loader do
  def prep_tail(project_name, bucket_name, lines)
    Etna::Clients::Metis::TailResponse.new(project_name, bucket_name, lines)
  end

  def labors_config(scripts, bucket_name="pics", autolink=true)
    JSON.parse({
      project_name: "labors",
      config: {
        bucket_name: bucket_name,
        autolink: autolink,
        models: {
          victim: {
            scripts: scripts
          }
        }
      }
    }.to_json)
  end
  context 'files' do
    it 'matches tails to files' do
      tail = prep_tail('labors', 'pics', [
        {"type":"file","id":50,"parent_id":nil,"node_name":"ignore.deceased.png","updated_at":"2023-08-03 22:39:17 +0000","file_hash":"0cc175b9c0f1b6a831c399e269772661","archive_id":nil},
        {"type":"file","id":51,"parent_id":40,"node_name":"LABORS-LION-H2-C1.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"file","id":52,"parent_id":40,"node_name":"NO-ID.deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"file","id":53,"parent_id":40,"node_name":"nonmatch_deceased.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"parent","id":40,"parent_id":nil,"node_name":"victims","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil},
        {"type":"parent","id":41,"parent_id":40,"node_name":"archived","updated_at":"2024-02-19 22:39:17 +0000","file_hash":nil,"archive_id":nil},
      ])

      config = labors_config([
        {
          type: "file",
          folder_path: "victims",
          file_match: "*.deceased*.png",
          attribute_name: "photo_deceased"
        }
      ])

      rules={
          victim: '^LABORS-LION-H\d+-C\d+$'
      }

      update = Metis::Loader.new(config, rules).update_for(tail)

      expect(update['revisions']).to eq(
          victim: {
              'LABORS-LION-H2-C1' => {
                  'photo_deceased' => {
                      original_filename: 'LABORS-LION-H2-C1.deceased.png',
                      path: 'metis://labors/pics/victims/LABORS-LION-H2-C1.deceased.png'
                  }
              }
          }
      )
    end
  end

  context 'file collections' do
    it 'matches tails to file collection scripts' do
      tail = prep_tail('labors', 'pics', [
        {"type":"file","id":50,"parent_id":nil,"node_name":"ignore.png","updated_at":"2023-08-03 22:39:17 +0000","file_hash":"0cc175b9c0f1b6a831c399e269772661","archive_id":nil},
        {"type":"file","id":51,"parent_id":41,"node_name":"family_photo.1.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"file","id":52,"parent_id":41,"node_name":"family_photo.2.png","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":nil},
        {"type":"parent","id":40,"parent_id":nil,"node_name":"family","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil},
        {"type":"parent","id":41,"parent_id":40,"node_name":"LABORS-LION-H2-C1","updated_at":"2024-02-19 22:39:17 +0000","file_hash":nil,"archive_id":nil},
      ])

      config = labors_config([
        {
          "type": "file_collection",
          "folder_path": "family",
          "file_match": "{LABORS-LION-H2-C1,LABORS-LION-H2-C2}/family_photo.*.png",
          "attribute_name": "family_photos"
        }
      ])

      rules = {
        victim: '^LABORS-LION-H\d+-C\d+$'
      }

      update = Metis::Loader.new(config, rules).update_for(tail)

      expect(update['revisions']).to eq(
        victim: {
          'LABORS-LION-H2-C1' => {
            'family_photos' => [
              {
                original_filename: 'family_photo.1.png',
                path: 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.1.png'
              },
              {
                original_filename: 'family_photo.2.png',
                path: 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.2.png'
              }
            ]
          }
        }
      )
    end
  end

  context 'data frames' do

    let(:tail) do
      prep_tail('labors', 'pics', [
            {"type":"file","id":51,"parent_id":41,"node_name":"village-1.tsv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
            {"type":"file","id":52,"parent_id":41,"node_name":"village-2.csv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":nil},
            {"type":"file","id":53,"parent_id":42,"node_name":"village-3.csv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":nil},
            {"type":"parent","id":41,"parent_id":nil,"node_name":"villages","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil},
            {"type":"parent","id":42,"parent_id":nil,"node_name":"other-villages","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil}
        ])
    end

    let(:model_std) do
      Etna::Clients::Magma::Models.new({
        'victim' => {
         'template' => {
          'identifier' => 'name',
          'parent' => 'monster',
          'attributes' => {
            'monster' => { 'attribute_type' => 'parent' },
            'name' =>  { 'attribute_type' =>  'identifier' },
            'species' =>  { 'attribute_type' =>  'string' },
            'target_name' =>  { 'attribute_type' =>  'string' }
          }
         }
        },
        'monster' => {
         'template' => {
          'identifier' => 'name',
          'parent' => 'stub',
          'attributes' => {
            'name' =>  { 'attribute_type' =>  'identifier' },
            'victim' =>  { 'attribute_type' =>  'collection' }
          }
         }
        }
      })
    end

    let(:model_table) do
      Etna::Clients::Magma::Models.new({
        'victim' => {
          'template' => {
            'identifier' => 'id',
            'parent' => 'monster',
            'attributes' => {
              'monster' => { 'attribute_type' => 'parent' },
              'species' => { 'attribute_type' => 'string' },
              'target_name' => { 'attribute_type' => 'string' }
            }
          },
        },
        'monster' => {
         'template' => {
          'identifier' => 'name',
          'parent' => 'stub',
          'attributes' => {
            'name' =>  { 'attribute_type' =>  'identifier' },
            'victim' =>  { 'attribute_type' =>  'table' }
          }
         }
        }
      })
    end

    let(:rules) do
      {
        monster: '^LABORS-LION$',
        victim: '^LABORS-LION-H\d+-C\d+$'
      }
    end

    let(:metis) do
      @metis
    end

    before(:each) do
      stub_download_file(project: 'labors', file: "pics/villages/village-1.tsv", file_contents: "name\tSPECIES\ttarget_name\nLABORS-LION-H1-C1\t__\t\n")
      stub_download_file(project: 'labors', file: 'pics/villages/village-2.csv', file_contents: "name,SPECIES,target_name\nLABORS-LION-H2-C1,__,\n")
      stub_download_file(project: 'labors', file: 'pics/other-villages/village-3.csv', file_contents: "monster,SPECIES,target_name\nLABORS-LION,human,Opheltes\n")
      stub_metis_setup
      @metis = Etna::Clients::Metis.new(host: METIS_HOST, token: TEST_TOKEN)
    end

    it 'updates from a tsv' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.tsv",
          "format": "tsv",
          "blank_table": true,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])

      update = Metis::Loader.new(config, rules, {}, model_std).update_for(
        tail, metis
      )

      expect(update['revisions']).to eq(
        'victim' => {
          'LABORS-LION-H1-C1' => {
            'species': '__',
            #'target_name': ''
          }
        }
      )
    end

    it 'updates from a csv' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.csv",
          "format": "csv",
          "blank_table": false,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])

      update = Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)

      expect(update['revisions']).to eq(
        'victim' => {
          'LABORS-LION-H2-C1' => {
            'species': '__',
            #'target_name': ''
          }
        }
      )
    end

    it 'auto-detects format and updates' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.{t,c}sv",
          "format": "auto-detect",
          "blank_table": false,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])

      update = Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)

      expect(update['revisions']).to eq(
        'victim' => {
          'LABORS-LION-H1-C1' => {
            'species': '__',
            #'target_name': ''
          },
          'LABORS-LION-H2-C1' => {
            'species': '__',
            #'target_name': ''
          }
        }
      )
    end

    it 'complains about missing columns in column_map' do
      config = labors_config(
        [
          {
            "type": "data_frame",
            "folder_path": "villages",
            "file_match": "*.tsv",
            "format": "auto-detect",
            "blank_table": false,
            "column_map": {
              "name": "name",
              "species": "species",
              "target_name": "target_name"
            }
          }
        ]
      )
      expect {
        Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)
      }.to(
        raise_error(Metis::Loader::Error, /missing column.*species/)
      )
    end

    it 'complains if there is no identifier attribute in the column map' do
      config = labors_config(
        [
          {
            "type": "data_frame",
            "folder_path": "villages",
            "file_match": "*.tsv",
            "format": "auto-detect",
            "blank_table": false,
            "column_map": {
              "nomen": "name",
              "species": "SPECIES",
              "target_name": "target_name"
            }
          }
        ]
      )
      expect {
        Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)
      }.to(
        raise_error(Metis::Loader::Error, /Identifier attribute is missing/)
      )
    end

    it 'complains if non-existent attributes are in the column_map' do
      config = labors_config(
        [
          {
            "type": "data_frame",
            "folder_path": "villages",
            "file_match": "*.tsv",
            "format": "auto-detect",
            "blank_table": false,
            "column_map": {
              "name": "name",
              "spoocies": "SPECIES",
              "target_name": "target_name"
            }
          }
        ]
      )

      # Now adjusting the model template to where column_map targets non-existent attributes
      expect{Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)}.to(
        raise_error(Metis::Loader::Error, /attribute\(s\) that don't exist: spoocies/)
      )
    end

    context 'table attributes' do
      it 'blanks a table attribute' do
          config = labors_config(
                          [
                              {
                                  "type": "data_frame",
                                  "folder_path": "other-villages",
                                  "file_match": "*.csv",
                                  "format": "csv",
                                  "blank_table": true,
                                  "column_map": {
                                      "monster": "monster",
                                      "species": "SPECIES",
                                      "target_name": "target_name"
                                  }
                              }
                          ]
          )
        update = Metis::Loader.new(config, rules, {}, model_table).update_for(tail, metis)

        expect(update['revisions']).to eq(
          'victim' => {
            '::temp-id-0' => {
              'monster': 'LABORS-LION',
              'species': 'human',
              'target_name': 'Opheltes'
            }
          },
          'monster' => {
            'LABORS-LION' => {'victim': ['::temp-id-0']}
          }
        )
      end

      it 'complains if a table has no parent attribute in column map' do
          config = labors_config(
            [
              {
                "type": "data_frame",
                "folder_path": "villages",
                "file_match": "*.tsv",
                "format": "auto-detect",
                "blank_table": false,
                "column_map": {
                  "munster": "monster",
                  "species": "SPECIES",
                  "target_name": "target_name"
                }
              }
            ]
          )
          expect{Metis::Loader.new(config, rules, {}, model_table).update_for(tail, metis)}.to(
            raise_error(Metis::Loader::Error, /Parent attribute is missing/)
          )
      end
    end

    it 'ignores hole values' do
      config = labors_config( [
          {
              "type": "data_frame",
              "folder_path": "villages",
              "file_match": "*.tsv",
              "format": "tsv",
              "blank_table": true,
              "column_map": {
                  "name": "name",
                  "species": "SPECIES",
                  "target_name": "target_name"
              },
              "values_to_ignore": "__,"
          }
        ]
      )
      update = Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)
      expect(update['revisions']).to eq(
        'victim' => {
          'LABORS-LION-H1-C1' => {}
        }
      )
    end

    it 'passes through dry_run and autolink' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.tsv",
          "format": "tsv",
          "blank_table": true,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])
      update = Metis::Loader.new(config, rules, {commit: false}, model_std).update_for(tail, metis)
      expect(update['dry_run']).to be_truthy
      expect(update['autolink']).to be_truthy
    end

    it 'complains about format' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.csv",
          "format": "tsv",
          "blank_table": false,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])
      expect {
        Metis::Loader.new(config, rules, {}, model_std).update_for(tail, metis)
      }.to(
        raise_error(Metis::Loader::Error, /Check the 'format' configuration/)
      )
    end
  end

end
