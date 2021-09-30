describe Polyphemus::Ipi::IpiWatchFoldersEtl do
  let(:etl) { Polyphemus::Ipi::IpiWatchFoldersEtl.new }
  let(:files) { [] }
  let(:project_name) { "ipi" }
  let(:folder_cursor) {
                        Polyphemus::MetisFolderEtlCursor.new(
                          job_name: "test",
                          project_name: project_name,
                          bucket_name: bucket_name,
                        )
                      }
  let(:cursor) { folder_cursor }
  let(:helper) { IpiHelper.new("lib/etls/renaming/projects/test_renames.json") }

  before(:each) do
    allow(IpiHelper).to receive(:new).and_return(helper)

    stub_metis_setup
    copy_renaming_project

    stub_magma_update_json
    stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")

    stub_request(:post, /#{METIS_HOST}\/#{project_name}\/find\/#{bucket_name}/)
      .to_return({
        status: 200,
        headers: {
          'Content-Type' => 'application/json'
        },
        body: {
                files: files,
                folders: [],
              }.to_json
      })
  end

  describe 'process_bulk_rna_seq_results' do
    let(:bucket_name) { 'data' }

    describe "processes files in found folders" do
      let(:files) do
        [{
          file_name: "rnaseq_table.tsv",
          download_url: "/#{project_name}/download/rnaseq_table.tsv",
          bucket_name: bucket_name,
          project_name: project_name,
          folder_id: 1,
        }]
      end

      describe "for attributes" do
        it 'works' do
          stub_download_file(
            project: project_name,
            file_contents: ::File.read("./spec/fixtures/ipi_rna_seq_results.tsv"),
          )

          etl.process(cursor, [
            create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
            create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
          ])

          # Make sure rna_seq records are updated, once per rna_seq record. Only when file is found.
          #   Three valid records in the given response, so three updates.
          expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(3)
        end
      end

      describe "for matrices" do
        let(:files) do
          [{
            file_name: "gene_counts_table.tsv",
            download_url: "/#{project_name}/download/gene_counts_table.tsv",
            bucket_name: bucket_name,
            project_name: project_name,
            folder_id: 1,
          }]
        end

        it 'works' do
          stub_download_file(
            project: project_name,
            file_contents: ::File.read("./spec/fixtures/ipi_gene_counts.tsv"),
          )

          etl.process(cursor, [
            create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
            create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
          ])

          # Make sure rna_seq records are updated, once per rna_seq record. Only when file is found.
          #   two records in the given response, so two updates.
          expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(2)
        end
      end
    end

    it "for invalid NASH / NAFLD samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl.process(cursor, [
        create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
        create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)
    end
  end

  describe 'link_rna_seq_raw_fastq_files' do
    let(:bucket_name) { "integral_data" }
    let(:files) do
      [
        {
          file_name: "something.fastq.gz",
          folder_id: 1,
          bucket_name: bucket_name,
          project_name: project_name,
        }
      ]
    end

    describe "create Polyphemus::WatchFile records" do
      it "for invalid NASH / NAFLD samples" do
        expect(Polyphemus::WatchFolder.count).to eq(0)

        etl.process(cursor, [
          create_metis_folder("IPIADR001.NASH1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NASH1.rna.live", id: 1),
          create_metis_folder("IPIADR001.NAFLD1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NAFLD1.rna.live", id: 2),
        ])

        expect(Polyphemus::WatchFolder.count).to eq(2)
      end

      it "for incorrectly named samples" do
        expect(Polyphemus::WatchFolder.count).to eq(0)

        etl.process(cursor, [
          create_metis_folder("WRONG001.T1.rna.tumor", "some_folder/BulkRNASeq/WRONG001.T1.rna.tumor", id: 1),
        ])

        expect(Polyphemus::WatchFolder.count).to eq(1)
      end
    end

    it "links files in found folders" do
      etl.process(cursor, [
        create_metis_folder("PATIENT001.T1.comp", "some_folder/BulkRNASeq/PATIENT001.T1.comp", id: 1),
        create_metis_folder("PATIENT001.N1.comp", "some_folder/BulkRNASeq/PATIENT001.N1.comp", id: 2),
        create_metis_folder("PATIENT002.T1.comp", "some_folder/BulkRNASeq/PATIENT002.T1.comp", id: 3),
      ])

      # Make sure rna_seq records are updated. Once per folder with files.
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
    end

    describe "create Polyphemus::WatchFile records" do
      it "for invalid NASH / NAFLD samples" do
        expect(Polyphemus::WatchFolder.count).to eq(0)

        etl.process(cursor, [
          create_metis_folder("IPIADR001.NASH1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NASH1.rna.live", id: 1),
          create_metis_folder("IPIADR001.NAFLD1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NAFLD1.rna.live", id: 2),
        ])

        expect(Polyphemus::WatchFolder.count).to eq(2)
      end

      it "for incorrectly named samples" do
        expect(Polyphemus::WatchFolder.count).to eq(0)

        etl.process(cursor, [
          create_metis_folder("WRONG001.T1.rna.tumor", "some_folder/BulkRNASeq/WRONG001.T1.rna.tumor", id: 1),
        ])

        expect(Polyphemus::WatchFolder.count).to eq(1)
      end
    end

    it "links files in found folders" do
      etl.process(cursor, [
        create_metis_folder("PATIENT001.T1.comp", "some_folder/BulkRNASeq/PATIENT001.T1.comp", id: 1),
        create_metis_folder("PATIENT001.N1.comp", "some_folder/BulkRNASeq/PATIENT001.N1.comp", id: 2),
        create_metis_folder("PATIENT002.T1.comp", "some_folder/BulkRNASeq/PATIENT002.T1.comp", id: 3),
      ])

      # Make sure rna_seq records are updated. Once per folder with files.
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
    end
  end

  describe 'link_processed_rna_seq_files' do
    let(:bucket_name) { "data" }

    let(:files) do
      [
        {
          file_name: "something.fastq.gz",
          folder_id: 1,
          bucket_name: bucket_name,
          project_name: project_name,
        }
      ]
    end

    describe 'magma record creation' do
      def folder(folder_name, folder_path, updated_at = Time.now)
        @folder_id ||= 0
        Etna::Clients::Metis::Folder.new({
          folder_name: folder_name,
          folder_path: folder_path,
          updated_at: updated_at,
          id: (@folder_id += 1)
        })
      end

      it "for all rna_seq" do
        etl.process(cursor, [
          folder("IPIADR001.N1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output/IPIADR001.N1.rna.live"),
          folder("IPIADR001.T1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output2/IPIADR001.T1.rna.live"),
          folder("IPIBLAD001.T1.rna.live", "bulkRNASeq/plate2_rnaseq_new/output/IPIBLAD001.T1.rna.live"),
        ])

        # Make sure plates are created
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq_plate": {
                "Plate1": {
                  "project": "UCSF Immunoprofiler",
                },
                "Plate2": {
                  "project": "UCSF Immunoprofiler",
                },
              },
            },
          }))

        # Make sure rna_seq records are created
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "IPIADR001.N1.rna.live": {
                  "rna_seq_plate": "Plate1",
                  "sample": "IPIADR001.N1",
                },
              },
            },
          }))
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "IPIBLAD001.T1.rna.live": {
                  "rna_seq_plate": "Plate2",
                  "sample": "IPIBLAD001.T1",
                },
              },
            },
          }))
      end

      it "does not create NASH / NAFLD samples" do
        etl.process(cursor, [
          folder("IPIADR001.NASH1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output/IPIADR001.NASH1.rna.live"),
          folder("IPIADR001.NAFLD1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output/IPIADR001.NAFLD1.rna.live"),
        ])

        # Plates created anyways ... no plate is purely ignored samples,
        #   so this is okay.
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq_plate": {
                "Plate1": {
                  "project": "UCSF Immunoprofiler",
                },
              },
            },
          }))

        # Make sure NO rna_seq records are created
        expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "IPIADR001.NASH1.rna.live": {
                  "rna_seq_plate": "Plate1",
                  "sample": "IPIADR001.NASH1",
                },
              },
            },
          }))
        expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "IPIADR001.NAFLD1.rna.live": {
                  "rna_seq_plate": "Plate1",
                  "sample": "IPIADR001.NAFLD1",
                },
              },
            },
          }))
      end

      it "for control" do
        etl.process(cursor, [
          folder("CONTROL_jurkat.plate1", "bulkRNASeq/plate1_rnaseq_new/output/CONTROL_jurkat.plate1"),
          folder("CONTROL_uhr.plate2", "bulkRNASeq/plate2_rnaseq_new/output/CONTROL_uhr.plate2"),
        ])

        # Make sure plates are created
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq_plate": {
                "Plate1": {
                  "project": "UCSF Immunoprofiler",
                },
                "Plate2": {
                  "project": "UCSF Immunoprofiler",
                },
              },
            },
          }))

        # Make sure Control record names work with validation and not attached to sample
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "Control_Jurkat.Plate1": {
                  "rna_seq_plate": "Plate1",
                },
              },
            },
          }))
        expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
          .with(body: hash_including({
            "revisions": {
              "rna_seq": {
                "Control_UHR.Plate2": {
                  "rna_seq_plate": "Plate2",
                },
              },
            },
          }))
      end
    end

    describe "create Polyphemus::WatchFile records" do
      it "for invalid NASH / NAFLD samples" do
        expect do
          etl.process(cursor, [
            create_metis_folder("IPIADR001.NASH1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output/IPIADR001.NASH1.rna.live", id: 1),
            create_metis_folder("IPIADR001.NAFLD1.rna.live", "bulkRNASeq/plate1_rnaseq_new/output/IPIADR001.NAFLD1.rna.live", id: 2),
          ])
        end.to change { Polyphemus::WatchFolder.count }.to(2)
      end

      it "for incorrectly named samples" do
        expect(Polyphemus::WatchFolder.count).to eq(0)

        etl.process(cursor, [
          create_metis_folder("WRONG001.T1.rna.tumor", "bulkRNASeq/plate1_rnaseq_new/output/WRONG001.T1.rna.tumor", id: 1),
        ])

        expect(Polyphemus::WatchFolder.count).to eq(1)
      end
    end

    it "links files in found folders" do
      etl.process(cursor, [
        create_metis_folder("PATIENT001.T1.comp", "bulkRNASeq/plate1_rnaseq_new/output/PATIENT001.T1.comp", id: 1),
        create_metis_folder("PATIENT001.N1.comp", "bulkRNASeq/plate1_rnaseq_new/output/PATIENT001.N1.comp", id: 2),
        create_metis_folder("PATIENT002.T1.comp", "bulkRNASeq/plate2_rnaseq_new/output/PATIENT002.T1.comp", id: 3),
      ])
      # Make sure rna_seq records are created and updated.
      # Once per process per folder when files found.
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(2)
    end
  end
end