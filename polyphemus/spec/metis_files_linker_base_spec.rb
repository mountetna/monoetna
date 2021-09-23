require_relative "./spec_helper"
require_relative '../lib/etls/metis_files_linker_base'

describe Polyphemus::MetisFilesLinkerBase do
  class MockMagmaClient < Etna::Clients::Magma
    def initialize()
      # does not call super to avoid validations made there.  This client is only meant to be
      # mocked and not invoked.
    end
  end

  let(:metis_folders) { {} }
  def make_file(
    file_path,
    file_name: File.basename(file_path),
    add_watch_folder: true
  )
    dir = File.dirname(file_path)

    unless metis_folders.include?(dir)
      folder_id = metis_folders[dir] = next_i
      create(
        :watch_folder,
        project_name: project_name,
        bucket_name: bucket_name,
        updated_at: "2021-01-01 00:00:00",
        folder_path: dir,
        watch_type: "link_files",
        metis_id: folder_id
      ) if add_watch_folder
    else
      folder_id = metis_folders[dir]
    end

    files << Etna::Clients::Metis::File.new({
      project_name: project_name,
      bucket_name: bucket_name,
      file_path: file_path,
      file_name: file_name,
      folder_id: folder_id,
    })
  end

  def next_i
    @i += 1
  end

  def attribute_regex
    {
      "one": /one.*$/,
      "two": /two.*$/,
    }
  end

  def make_record(attrs)
    record_name = "record-#{next_i}"
    existing_records[record_name] = attrs
    return [record_name, attrs]
  end

  let(:model_name) { 'somemodel' }
  let(:magma_client) { MockMagmaClient.new }
  let(:project_name) { "BLARG" }
  let(:bucket_name) { "bucketboy" }
  let(:existing_records) { }
  let(:record_name_regex) do
    /.*\/(?<record_name>.*)\/.*$/
  end
  let(:files) { [ ] }
  let(:existing_records) { { } }

  let(:files_by_record_name) do
    linker.organize_metis_files_by_magma_record(
      metis_files: files,
      magma_record_names: existing_records.keys,
      path_regex: record_name_regex,
    )
  end

  let(:linker) do
    Polyphemus::MetisFilesLinkerBase.new(project_name: project_name, bucket_name: bucket_name)
  end

  before(:each) do
    @i = 0

    allow(linker).to receive(:magma_client).and_return(magma_client)
    allow(magma_client).to receive(:update_json)
    allow(magma_client).to receive(:retrieve).and_return(
      Etna::Clients::Magma::RetrievalResponse.new({
        'models' => {
          model_name => {
            'documents' => existing_records,
            'template' => {
              'attributes' => {
                'one' => {
                  'attribute_type' => 'file',
                },
                'two' => {
                  'attribute_type' => 'file_collection',
                },
              }
            }
          }
        }
      })
    )
  end


  describe 'for existing records' do
    let(:with_backing_folder_watch) { true }
    describe 'but without a file watch backing it' do
      it 'runs without exception, processing nothing' do
        @record_one_name, attrs = make_record({})
        make_file("processed/#{@record_one_name}/one-a", add_watch_folder: false)

        expect(files_by_record_name).to be_empty
        expect(linker.link(
          model_name: model_name,
          files_by_record_name: files_by_record_name,
          attribute_regex: attribute_regex
        )).to eql(false)
      end
    end

    it 'can merge in a few files into a non empty collection attribute' do
      @record_one_name, attrs = make_record({
        two: [
          {
            original_filename: "two-a",
            path: "metis://#{project_name}/#{bucket_name}/processed/record-1/two-a"
          },
          {
            original_filename: "two-b",
            path: "metis://#{project_name}/#{bucket_name}/processed/record-1/two-b"
          },
        ]
      })
      make_file("processed/#{@record_one_name}/two-a")
      # two-b is no longer included, but should not be dropped
      make_file("processed/#{@record_one_name}/two-c")

      expect(linker.link(
        model_name: model_name,
        files_by_record_name: files_by_record_name,
        attribute_regex: attribute_regex
      )).to eql({
        model_name => {
          @record_one_name => {
            two: [
              {
                original_filename: "two-a",
                path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/two-a"
              },
              {
                original_filename: "two-c",
                path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/two-c"
              },
              {
                original_filename: "two-b",
                path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/two-b"
              },
            ]
          }
        }
      })
    end

    it 'can link in a few files into an empty collection attribute' do
      @record_one_name, attrs = make_record({})
      make_file("processed/#{@record_one_name}/two-a")
      make_file("processed/#{@record_one_name}/two-b")

      expect(linker.link(
        model_name: model_name,
        files_by_record_name: files_by_record_name,
        attribute_regex: attribute_regex
      )).to eql({
        model_name => {
          @record_one_name => {
            two: [
              {
                original_filename: "two-a",
                path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/two-a"
              },
              {
                original_filename: "two-b",
                path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/two-b"
              },
            ]
          }
        }
      })
    end

    it 'can link in a single file attribute' do
      @record_one_name, attrs = make_record({})
      make_file("processed/#{@record_one_name}/one-a")

      expect(files_by_record_name).to_not be_empty

      expect(linker.link(
        model_name: model_name,
        files_by_record_name: files_by_record_name,
        attribute_regex: attribute_regex
      )).to eql({
        model_name => {
          @record_one_name => {
            one: {
              original_filename: "one-a",
              path: "metis://#{project_name}/#{bucket_name}/processed/#{@record_one_name}/one-a"
            }
          }
        }
      })
    end
  end
end