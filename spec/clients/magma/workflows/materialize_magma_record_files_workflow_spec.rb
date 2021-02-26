require 'webmock/rspec'
require 'fileutils'
require 'json'

describe Etna::Clients::Magma::MaterializeDataWorkflow do
  describe "e2e" do
    MATERIALIZE_ATTRIBUTES = {
        'patient' => [
            'ice_time',
            'phyician',
            'sample',
        ],
        'sample' => [
            'population',
            'rna_seq',
            'dc_file',
        ],
        'rna_seq' => [
            'tube_name',
            'rna_seq_plate',
            'gene_counts',
        ],
        'rna_seq_plate' => [
            'submission_date',
            'rna_seq',
        ],
    }
    xit 'can materialize metadata and files into a file system class' do
      project_name = 'ipi'
      configure_etna_yml

      VCR.use_cassette('materialize_data_workflow-base_case.e2e') do
        magma_client = Etna::Clients::Magma.new(
            host: 'https://magma.development.local',
            token: ENV['TOKEN'] || TEST_TOKEN,
            ignore_ssl: true,
        )

        metis_client = Etna::Clients::Metis.new(
            host: 'https://metis.development.local',
            token: ENV['TOKEN'] || TEST_TOKEN,
            ignore_ssl: true,
        )

        allow(metis_client).to receive(:file_metadata).and_return({ size: 100, etag: 'the-returned-hash' })

        workflow = Etna::Clients::Magma::MaterializeDataWorkflow.new(
            metis_client: metis_client,
            magma_client: magma_client,
            project_name: project_name,
            model_name: 'patient',
            model_filters: {
                'patient' => 'ice_time>1',
                'rna_seq' => 'eisenberg_score>8',
            },
            model_attributes_mask: MATERIALIZE_ATTRIBUTES,
            filesystem: Etna::Filesystem.new,
            stub_files: true,
        )

        workflow.materialize_all do |dir|
          entries_by_dir = ::Dir.glob("**/*", base: dir).inject({}) do |acc, n|
            acc.tap do
              dirname = ::File.dirname(n)
              acc[dirname] ||= 0
              acc[dirname] += 1
            end
          end

          expect(entries_by_dir).to eq({
              "." => 6,
              "bin" => 1,
              "patient" => 160,
              "rna_seq" => 258,
              "rna_seq_plate" => 7,
              "sample" => 287,
          })
        end
      end
    end
  end
end
