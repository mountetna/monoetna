describe Polyphemus::CascadeMvirPatientWaiverToRestricted do
  let(:command) { Polyphemus.instance.subcommands['cascade_mvir_patient_waiver_to_restricted'] }

  context 'for restricting patients and pools' do
    let(:patients_to_restrict) {
      {
        'Sally': { name: 'Sally', restricted: false, consent: 'Full Study' },
        'Dan': { name: 'Dan', restricted: false, consent: 'Initial Waiver' },
        'Mike': { name: 'Mike', restricted: false, consent: 'No Samples/Data' }
      }
    }

    before(:each) do
      @all_updates = []

      WebMock.disable_net_connect!
      stub_magma_setup(patients_to_restrict)
      stub_magma_restricted_pools('cytof', ['pool-a', 'pool-c'])
      stub_magma_all_pools('cytof', ['pool-a', 'pool-c'])
      stub_magma_restricted_pools('sc_rna_seq', ['pool-e', 'pool-f'])
      stub_magma_all_pools('sc_rna_seq', ['pool-e', 'pool-f'])
      stub_metis_setup
    end

    after(:each) do
      WebMock.enable_net_connect!
    end


    it 'updates the restricted states correctly for pools and the patients' do
      stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_folder({bucket: RELEASE_BUCKET})
      stub_bucket_find(
        bucket: RELEASE_BUCKET,
        response_body: {
          folders: []
        }
      )
      stub_bucket_find(
        bucket: RESTRICT_BUCKET,
        response_body: {
          folders: []
        }
      )

      command.execute

      expect(@all_updates.inject({}) do |acc, n|
        n.keys.each do |k|
          (acc[k] ||= {}).update(n[k])
        end
        acc
      end).to eq({
        "cytof_pool" => {
          "pool-a" => {"restricted" => "true"},
          "pool-c" => {"restricted" => "true"}
        },
          "sc_rna_seq_pool" => {
              "pool-e" => {"restricted" => "true"},
              "pool-f" => {"restricted" => "true"},
          },
        "patient" => {
            "Mike" => {"restricted" => "true"},
            "Dan" => {"restricted" => "true"},
        }
      })
    end

    it 'correctly makes Metis calls to restrict patient and pool folders' do
      patients_to_restrict[:Dan][:consent] = 'No Samples/Data'
      stub_magma_setup(patients_to_restrict)

      stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_folder({bucket: RELEASE_BUCKET})
      stub_bucket_find(
        bucket: RELEASE_BUCKET,
        response_body: {
          folders: [
            create_metis_folder("Dan-D0-ASSAY1", "assay/processed/Dan-D0-ASSAY1").raw,
            create_metis_folder("Dan-D12-ASSAY1", "assay/processed/Dan-D12-ASSAY1").raw,
            create_metis_folder("Dan-D0-ASSAY3", "assay/raw/Dan-D0-ASSAY3").raw,
            create_metis_folder("Dan-D12-ASSAY3", "assay/raw/Dan-D12-ASSAY3").raw,
          ]
        },
        response_body_2: {
          folders: [
            create_metis_folder("Mike-D0-ASSAY3", "assay/processed/Mike-D0-ASSAY3").raw,
            create_metis_folder("Mike-D12-ASSAY3", "assay/processed/Mike-D12-ASSAY3").raw,
            create_metis_folder("Mike-D0-ASSAY", "assay/raw/Mike-D0-ASSAY").raw,
            create_metis_folder("Mike-D12-ASSAY4", "assay/raw/Mike-D12-ASSAY4").raw
          ]
        },
        response_body_3: {
          folders: [
            create_metis_folder("pool-a", "assay/processed/pool-a").raw,
            create_metis_folder("pool-a", "assay/raw/pool-a").raw,
          ]
        }
      )
      stub_bucket_find(
        bucket: RESTRICT_BUCKET,
        response_body: {
          folders: []
        }
      )

      command.execute

      # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
      # This should ignore the "summary" sub-folder that is under a patient
      # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
      # This should ignore the "summary" sub-folder that is under the pool
      # pool-c not in the release folder fixture, so it shouldn't have requests
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).times(10)
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).times(10)
    end

    it 'continues working if a single patient throws an error' do
      stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_folder_with_error({bucket: RELEASE_BUCKET})
      stub_bucket_find(
        bucket: RELEASE_BUCKET,
        response_body: {
          folders: [
            create_metis_folder("Dan-D0-ASSAY1", "assay/processed/Dan-D0-ASSAY1").raw,
            create_metis_folder("Dan-D12-ASSAY1", "assay/processed/Dan-D12-ASSAY1").raw,
            create_metis_folder("Dan-D0-ASSAY3", "assay/raw/Dan-D0-ASSAY3").raw,
            create_metis_folder("Dan-D12-ASSAY3", "assay/raw/Dan-D12-ASSAY3").raw,
          ]
        },
        response_body_2: {
          folders: [
            create_metis_folder("Mike-D0-ASSAY3", "assay/processed/Mike-D0-ASSAY3").raw,
            create_metis_folder("Mike-D12-ASSAY3", "assay/processed/Mike-D12-ASSAY3").raw,
            create_metis_folder("Mike-D0-ASSAY", "assay/raw/Mike-D0-ASSAY").raw,
            create_metis_folder("Mike-D12-ASSAY4", "assay/raw/Mike-D12-ASSAY4").raw
          ]
        },
        response_body_3: {
          folders: [
            create_metis_folder("pool-a", "assay/processed/pool-a").raw,
            create_metis_folder("pool-a", "assay/raw/pool-a").raw,
          ]
        }
      )
      stub_bucket_find(
        bucket: RESTRICT_BUCKET,
        response_body: {
          folders: []
        }
      )

      command.execute

      # One patient skips, so only 7 requests per patient (fewer than above test)
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).times(7)
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).times(7)
    end
  end

  context 'when releasing patients and pools' do
    let(:patients_to_release) {
      {
        'Danielle': { name: 'Danielle', restricted: true, consent: 'Full Study' },
        'Dan': { name: 'Dan', restricted: false, consent: 'Full Study' }
      }
    }

    before(:each) do
      @all_updates = []

      WebMock.disable_net_connect!

      stub_magma_setup(patients_to_release)
      stub_magma_restricted_pools('cytof', [])
      stub_magma_restricted_pools('sc_rna_seq', [])
      stub_magma_all_pools('cytof', ['pool-b', 'pool-d'])
      stub_magma_all_pools('sc_rna_seq', ['pool-g', 'pool-j'])
      stub_metis_setup
    end

    it 'updates the restricted states correctly for pools and the patients' do
      stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})
      stub_create_folder({bucket: RELEASE_BUCKET})
      stub_rename_folder({bucket: RESTRICT_BUCKET})
      stub_bucket_find(
        bucket: RELEASE_BUCKET,
        response_body: {
          folders: []
        }
      )
      stub_bucket_find(
        bucket: RESTRICT_BUCKET,
        response_body: {
          folders: []
        }
      )

      command.execute

      expect(@all_updates.inject({}) do |acc, n|
        n.keys.each do |k|
          (acc[k] ||= {}).update(n[k])
        end
        acc
      end).to eq({
        "cytof_pool" => {
          "pool-b" => {"restricted" => "false"},
          "pool-d" => {"restricted" => "false"},
        },
        "sc_rna_seq_pool" => {
            "pool-g" => {"restricted" => "false"},
            "pool-j" => {"restricted" => "false"},
        },
        "patient" => {
            "Danielle" => {"restricted" => "false"},
        }
      })
    end

    it 'correctly makes Metis calls to release patient and pool folders' do
      stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})
      stub_create_folder({bucket: RELEASE_BUCKET})
      stub_rename_folder({bucket: RESTRICT_BUCKET})
      stub_bucket_find(
        bucket: RESTRICT_BUCKET,
        response_body: {
          folders: [
            create_metis_folder("Danielle-D7-ASSAY1", "assay/processed/Danielle-D7-ASSAY1").raw,
            create_metis_folder("Danielle-D14-ASSAY1", "assay/processed/Danielle-D14-ASSAY1").raw,
            create_metis_folder("Danielle-D7", "assay/raw/Danielle-D7").raw,
            create_metis_folder("Danielle-D14-ASSAY2", "assay/raw/Danielle-D14-ASSAY2").raw,
          ]
        },
        response_body_2: {
          folders: [
            create_metis_folder("pool-c", "assay/processed/pool-c").raw,
            create_metis_folder("pool-c", "assay/raw/pool-c").raw,
          ]
        }
      )
      stub_bucket_find(
        bucket: RELEASE_BUCKET,
        response_body: {
          folders: []
        }
      )

      command.execute

      # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
      # This should ignore the "summary" sub-folder that is under a patient
      # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
      # This should ignore the "summary" sub-folder that is under the pool
      # pool-d isn't in the restricted folder fixture, so it shouldn't have requests
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}/).times(6)
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/).times(6)

    end
  end
end
