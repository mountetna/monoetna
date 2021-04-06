require 'fileutils'

describe Vulcan::Storage do
  let(:storage) { Vulcan::Storage.new }
  before(:each) do
    FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)

    @input_file_a = Vulcan::Storage::StorageFile.new(project_name: 'project', cell_hash: 'abc', data_filename: 'r', logical_name: 'a', prefix: 'a')
    @input_file_b = Vulcan::Storage::StorageFile.new(project_name: 'project', cell_hash: 'abc', data_filename: 'r', logical_name: 'b', prefix: 'a')
    @input_files = [@input_file_a, @input_file_b]
    @output_filenames = ['a', 'b']
    @session_key = 'session_key'
    @script = 'script'
    @raw_hash = nil
    @project_name = 'project'
  end

  def build_target
    Vulcan::Storage::BuildTarget.new(project_name: @project_name, session_key: @session_key, input_files: @input_files, output_filenames: @output_filenames, script: @script)
  end

  def material_source
    Vulcan::Storage::MaterialSource.new(project_name: @project_name, session_key: @session_key, material_reference: {sha: "abc"})
  end

  describe '#with_build_transaction' do
    def test_with_build_transaction(buildable)
      expect do
        storage.with_run_cell_build(run_cell: storage.run_cell_for(buildable)) do |build_dir, output_files|
          output_files = output_files.map { |of| of.data_path(storage) }
          output_files.each { |of| expect(::File.exists?(of)).to eql(true) }
          output_files.each { |of| ::File.write(of, "test data") }
          output_files.each { |of| expect(::File.read(of)).to_not eql("") }

          raise "Oh no!"
        end
      end.to raise_error("Oh no!")

      built_paths = buildable.build_outputs.values.map { |sf| sf.data_path(storage) }
      built_paths.each { |of| expect(::File.exists?(of)).to eql(false) }

      storage.with_run_cell_build(run_cell: storage.run_cell_for(buildable)) do |build_dir, output_files|
        output_files = output_files.map { |of| of.data_path(storage) }

        output_files.each { |of| expect(::File.exists?(of)).to eql(true) }
        output_files.each { |of| expect(::File.read(of)).to eql("") }
        output_files.each_with_index { |of, i| ::File.write(of, "stuff #{i}") }

        storage.install_build_output(buildable: buildable, build_dir: build_dir)
      end

      built_paths.each { |of| expect(::File.exists?(of)).to eql(true) }
      built_paths.each_with_index { |of, i| expect(::File.read(of)).to eql("stuff #{i}") }
    end

    it 'works with BuildTarget' do
      test_with_build_transaction(build_target)
    end

    it 'works with MaterialSource' do
      test_with_build_transaction(material_source)
    end
  end

  describe 'BuildTarget' do

    describe '#cell_hash' do
      def cell_hash
        build_target.cell_hash
      end

      it 'is stable for input file ordering' do
        expect { @input_files.reverse! }.to_not change { cell_hash }
      end

      it 'is stable for input file prefix' do
        expect { @input_file_a.prefix += 'moo' }.to_not change { cell_hash }
      end

      it 'is stable for output file name ordering' do
        expect { @output_filenames.reverse! }.to_not change { cell_hash }
      end

      it 'changes by the input file hashes' do
        expect { @input_file_a.cell_hash += 'b' }.to change { cell_hash }
      end

      it 'changes when the input is a script or a raw_hash' do
        expect { @raw_hash = @script; @script = nil }.to change { cell_hash }
      end

      it 'changes when session_key changes' do
        expect { @session_key += 'session_key_more' }.to change { cell_hash }
      end
    end
  end

  describe '.hash_json_obj' do
    def hash(obj)
      Vulcan::Storage.hash_json_obj(obj)
    end

    it 'does not process non json objects' do
      expect { hash(Set.new) }.to raise_error
    end

    it 'does process storage file objects' do
      values = {project_name: 'abc', cell_hash: '123', data_filename: 'one', logical_name: 'logical_name'}
      expect(hash(Vulcan::Storage::StorageFile.new(**values))).to eql(hash(values))
    end

    it 'treats null entries for keys the same as undefined entries for keys' do
      expect(hash({a: 1, b: nil})).to eql(hash({a: 1}))
      expect(hash([nil])).to_not eql(hash([]))
    end

    it 'produces unique values per type' do
      expect(hash(1.0)).to eql(hash(1.0))
      expect(hash(1.0)).to_not eql(hash(1))
      expect(hash(1)).to eql(hash(1))
      expect(hash("1.0")).to_not eql(hash(1.0))
      expect(hash("1.0")).to eql(hash("1.0"))
      expect(hash(true)).to_not eql(hash(false))
      expect(hash(nil)).to_not eql(hash(false))
    end

    it 'gives unique objects stable identifiers' do
      expect(hash({
          a: 1
      })).to_not eql(hash({
          a: [1]
      }))


      expect(hash({
          a: [1],
          b: 2,
          c: {
              d: "3",
              e: true,
          }
      })).to eql(hash({
          b: 2,
          c: {
              e: true,
              d: "3",
          },
          a: [1],
      }))
    end
  end
end
