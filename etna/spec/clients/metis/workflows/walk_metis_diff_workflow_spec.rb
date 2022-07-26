
describe Etna::Clients::Metis::WalkMetisDiffWorkflow do
  def walker(metis_client:, project_name:, bucket_name:, root_dir:)
    Etna::Clients::Metis::WalkMetisWorkflow.new(
      metis_client: metis_client,
      project_name: project_name,
      bucket_name: bucket_name,
      root_dir: root_dir,
    )
  end

  describe 'unit tests' do
    def make_file(file_path:, updated_at:, file_hash:)
      Etna::Clients::Metis::File.new({
        file_path: file_path,
        updated_at: updated_at,
        file_hash: file_hash
      })
    end

    def make_folder(folder_path:, updated_at:)
      Etna::Clients::Metis::Folder.new({
        folder_path: folder_path,
        updated_at: updated_at,
      })
    end

    def path_of(file_or_folder)
      if file_or_folder.is_a?(Etna::Clients::Metis::Folder)
        file_or_folder.folder_path
      elsif file_or_folder.is_a?(Etna::Clients::Metis::File)
        file_or_folder.file_path
      end
    end

    def run_test(l, r)
      workflow = Etna::Clients::Metis::WalkMetisDiffWorkflow.new(
        left_walker: [[l, path_of(l)]],
        right_walker: [[r, path_of(r)]],
      )

      result = workflow.to_enum.to_a

      result.each do |row|
        expect(row[1].nil? || row[1] == l).to be_truthy
        expect(row[2].nil? || row[2] == r).to be_truthy
      end

      result.map(&:first)
    end

    it 'reports when one side is a folder' do
      expect(run_test(
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
      )).to eql([:right_is_folder])

      expect(run_test(
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
      )).to eql([:left_is_folder])
    end

    it 'reports missing files' do
      expect(run_test(
        make_folder(folder_path: '/c', updated_at: Time.at(5).to_s),
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
      )).to eql([:right_unique, :left_unique])

      expect(run_test(
        make_folder(folder_path: '/c', updated_at: Time.at(5).to_s),
        make_folder(folder_path: '/e', updated_at: Time.at(5).to_s),
      )).to eql([:left_unique, :right_unique])
    end

    it 'reports when two folders equal' do
      expect(run_test(
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
      )).to eql([:equal])

      expect(run_test(
        make_folder(folder_path: '/a', updated_at: Time.at(1).to_s),
        make_folder(folder_path: '/a', updated_at: Time.at(5).to_s),
      )).to eql([:left_older])

      expect(run_test(
        make_folder(folder_path: '/a', updated_at: Time.at(1).to_s),
        make_folder(folder_path: '/a', updated_at: nil),
      )).to eql([:unknown])
    end

    it 'reports when two files equal' do
      expect(run_test(
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
      )).to eql([:equal])

      expect(run_test(
        make_file(file_path: '/a', updated_at: Time.at(9).to_s, file_hash: 'a'),
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
      )).to eql([:equal])

      expect(run_test(
        make_file(file_path: '/a', updated_at: Time.at(9).to_s, file_hash: 'b'),
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: 'a'),
      )).to eql([:right_older])

      expect(run_test(
        make_file(file_path: '/a', updated_at: Time.at(9).to_s, file_hash: 'b'),
        make_file(file_path: '/a', updated_at: Time.at(5).to_s, file_hash: nil),
      )).to eql([:unknown])
    end
  end


  describe 'e2e' do
    it 'runs without error and completes' do
      VCR.use_cassette('walk_metis_diff_workflow.e2e') do
        metis_client = Etna::Clients::Metis.new(
          host: 'https://metis.ucsf.edu',
          token: ENV['TOKEN'] || TEST_TOKEN,
        )

        workflow = Etna::Clients::Metis::WalkMetisDiffWorkflow.new(
          left_walker: walker(
            metis_client: metis_client,
            project_name: 'mvir1',
            bucket_name: 'GNE_composite',
            root_dir: 'Upload/GNE_redacted_data/single_cell_ADT/raw',
          ),
          right_walker: walker(
            metis_client: metis_client,
            project_name: 'mvir1',
            bucket_name: 'GNE_redacted_data',
            root_dir: 'data/single_cell_ADT/raw',
          ),
        )
      end
    end
  end
end