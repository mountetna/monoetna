describe Metis::Loader do
  def prep_tail(project_name, bucket_name, lines)
    Etna::Clients::Metis::TailResponse.new(lines)
  end

  def labors_config(scripts, bucket_name="pics")
    {
      project_name: "labors",
      config: {
        bucket_name: bucket_name,
        autolink: true,
        models: {
          victim: {
            scripts: scripts
          }
        }
      }
    }
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

      # we have a config
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

      # given these we get an update
      update = Metis::Loader.new(config, rules).update_for(tail)

      require 'pry'
      binding.pry

      expect(update['revisions']).to eq(
          'victim': {
              'LABORS-LION-H2-C1': {
                  'photo_deceased': {
                      'original_filename': 'LABORS-LION-H2-C1.deceased.png',
                      'path': 'metis://labors/pics/victims/LABORS-LION-H2-C1.deceased.png'
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

      # we have a config
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

      # given these we get an update
      update = Metis::Loader.new(config, rules).update_for(tail)

      expect(update['revisions']).to eq(
        'victim': {
          'LABORS-LION-H2-C1': {
            'family_photos': [
              {
                'original_filename': 'family_photo.1.png',
                'path': 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.1.png'
              },
              {
                'original_filename': 'family_photo.2.png',
                'path': 'metis://labors/pics/family/LABORS-LION-H2-C1/family_photo.2.png'
              }
            ]
          }
        }
      )
    end
  end

  context 'data frames' do
    def Template(**r)
      r
    end
    def Attribute(r)
      r
    end
    let(:files) do
      {
            'village-1.tsv': 'name\tSPECIES\ttarget_name\nLABORS-LION-H2-C1\t__\t\n',
            'village-2.csv': 'name,SPECIES,target_name\nLABORS-LION-H2-C1,__,\n'
        }
    end

    let(:tail) do
      prep_tail('labors', 'pics', [
            {"type":"file","id":51,"parent_id":41,"node_name":"village-1.tsv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"8277e0910d750195b448797616e091ad","archive_id":nil},
            {"type":"file","id":52,"parent_id":41,"node_name":"village-2.csv","updated_at":"2023-11-11 22:39:17 +0000","file_hash":"9277e0910d750195b448797616e091ad","archive_id":nil},
            {"type":"parent","id":41,"parent_id":nil,"node_name":"villages","updated_at":"2023-08-03 22:39:17 +0000","file_hash":nil,"archive_id":nil}
        ])
    end
    let(:model_std) do
      {
            'victim': {
                'template': Template(
                    identifier:'name',
                    parent:'stub',
                    attributes:{
                        'name': Attribute({ 'attribute_type': 'identifier' }),
                        'species': Attribute({ 'attribute_type': 'string' }),
                        'target_name': Attribute({ 'attribute_type': 'string' })
                    }
                ),
                'isTable': false
            }
        }
    end
    let(:model_table) do
      {
            'victim': {
                'template': Template(
                    identifier:'id',
                    parent:'name',
                    attributes:{
                        'name': Attribute({ 'attribute_type': 'parent' }),
                        'species': Attribute({ 'attribute_type': 'string' }),
                        'target_name': Attribute({ 'attribute_type': 'string' })
                    }
                ),
                'isTable': true
            }
        }
    end

    it '' do
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

      update = Metis::Loader.new(config, @rules).update_for(
        tail, metis, model_std
      )

      expect(update['revisions']).to eq(
        'victim': {
          'LABORS-LION-H2-C1': {
            'species': '__',
            'target_name': ''
          }
        }
      )
    end

    it '' do
      config = labors_config([
        {
          "type": "data_frame",
          "folder_path": "villages",
          "file_match": "*.csv",
          "format": "auto-detect",
          "blank_table": false,
          "column_map": {
            "name": "name",
            "species": "SPECIES",
            "target_name": "target_name"
          }
        }
      ])
      update = Metis::Loader.new(config, @rules).update_for(tail, metis, model_std)
      expect(update['revisions']).to eq(
        'victim': {
          'LABORS-LION-H2-C1': {
            'species': '__',
            'target_name': ''
          }
        }
      )
    end

    it '' do
      config = labors_config(
        [
          {
            "type": "data_frame",
            "folder_path": "villages",
            "file_match": "*.tsv",
            "format": "auto-detect",
            "blank_table": true,
            "column_map": {
              "name": "name",
              "species": "species",
              "target_name": "target_name"
            }
          }
        ]
      )
      update = Metis::Loader.new(config, @rules).update_for(tail, metis, model_table)
      expect(update['revisions']).to eq(
        'name': {
          'LABORS-LION-H2-C1': {'victim': ['::temp-id-0']}
        },
        'victim': {
          '::temp-id-0': {
            'name': 'LABORS-LION-H2-C1',
            'species': '__',
            'target_name': ''
          }
        }
      )
    end

    it '' do
        config = labors_config(
                        [
                            {
                                "type": "data_frame",
                                "folder_path": "villages",
                                "file_match": "*.tsv",
                                "format": "csv",
                                "blank_table": true,
                                "column_map": {
                                    "name": "name",
                                    "species": "SPECIES",
                                    "target_name": "target_name"
                                }
                            }
                        ]
        )
      update = Metis::Loader.new(config, @rules).update_for(tail, metis, model_table)

      expect(update['revisions']).to eq(
        'victim': {
          '::temp-id-0': {
            'name': 'LABORS-LION-H2-C1',
            'species': '__',
            'target_name': ''
          }
        }
      )
    end

    it '' do
      config1hole = labors_config( [
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
      update = Metis::Loader.new(config1hole, @rules).update_for(tail, metis, model_std)
      expect(update['revisions']).to eq(
        'victim': {
          'LABORS-LION-H2-C1': {}
        }
      )
    end

    it '' do
      update = Metis::Loader.new(config1, @rules, params={"commit": false}).update_for(tail, metis, model_std)
      expect(update['dry_run']).to be_truthy
      expect(update['autolink']).to be_truthy
    end
  end

  it '' do
    expect(Metis::Loader.new(config3, @rules).update_for(tail, metis, model_std)).to(
      raise_error(Metis::LoaderError, /missing column.*species/)
    )

    expect(Metis::Loader.new(config3, @rules).update_for(tail, metis, model_table)).to(
      raise_error(Metis::LoaderError, /missing column.*species/)
    )

    expect(Metis::Loader.new(config4, @rules).update_for(tail, metis, model_std)).to(
      raise_error(Metis::LoaderError, /Check the 'format' configuration/)
    )

    expect(Metis::Loader.new(config4, @rules).update_for(tail, metis, model_table)).to(
      raise_error(Metis::LoaderError, /Check the 'format' configuration/)
    )

      # Now adjusting the model template to where the identifier would be missing from the column_map
      def models_id(isTable)
          return {
              'victim': {
                  'template': Template(
                      identifier='id',
                      attributes={
                          'name': Attribute({ 'attribute_type': 'identifier' }),
                          'species': Attribute({ 'attribute_type': 'string' }),
                          'target_name': Attribute({ 'attribute_type': 'string' })
                      }
                  ),
                  'isTable': isTable
              }
          }
      end
      expect(Metis::Loader.new(config1, @rules).update_for(tail, metis, models_id(false))).to(
        raise_error(Metis::LoaderError, /Identifier attribute is missing/)
      )
      expect(Metis::Loader.new(config1, @rules).update_for(tail, metis, models_id(true))).to(
        raise_error(Metis::LoaderError, /Parent attribute is missing/)
      )

      # Now adjusting the model template to where column_map targets non-existent attributes
      def models_missing(isTable)
          return {
              'victim': {
                  'template': Template(
                      identifier='name',
                      attributes={
                          'name': Attribute({ 'attribute_type': 'identifier' }),
                          'target_name': Attribute({ 'attribute_type': 'string' })
                      }
                  ),
                  'isTable': false
              }
          }
      end
      expect(Metis::Loader.new(config1, @rules).update_for(tail, metis, models_missing(false))).to(
        raise_error(Metis::LoaderError, /attribute\(s\) that don't exist: species/)
      )
      expect(Metis::Loader.new(config1, @rules).update_for(tail, metis, models_missing(true))).to(
        raise_error(Metis::LoaderError, /attribute\(s\) that don't exist: species/)
      )
  end
end
