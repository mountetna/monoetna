require_relative '../lib/folder_path_calculator'

describe Metis::FolderPathCalculator do

  before(:each) do
    @bucket = default_bucket('athena')

    @public_folder = create_folder('athena', 'public', bucket: @bucket)
    stubs.create_folder('athena', 'files', 'public')

    @child_folder = create_folder('athena', 'child', bucket: @bucket, folder: @public_folder)
    stubs.create_folder('athena', 'files', 'public/child')

    @grandchild_folder = create_folder('athena', 'grandchild', bucket: @bucket, folder: @child_folder)
    stubs.create_folder('athena', 'files', 'public/child/grandchild')
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'returns folder name for root folder' do
    calc = Metis::FolderPathCalculator.new(bucket: @bucket)
    expect(calc.get_folder_path(@public_folder.id)).to eq('public')
  end

  it 'correctly calculates path for nested folders' do
    calc = Metis::FolderPathCalculator.new(bucket: @bucket)
    expect(calc.get_folder_path(@grandchild_folder.id)).to eq('public/child/grandchild')
  end
end