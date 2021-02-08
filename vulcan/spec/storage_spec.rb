require 'fileutils'

describe Vulcan::Storage do
  describe '#with_build_transaction' do
    let(:storage) { Vulcan::Storage.new }
    before(:each) do
      FileUtils.rm_rf(storage.data_root)
    end

    it 'builds within a transaction' do
      expect do
        storage.with_build_transaction(project_name: 'project', ch: 'ch', output_filenames: ['a', 'b']) do |output_files|
          output_files = output_files.map { |of| of.data_path(storage) }
          expect(::File.exists?(output_files[0])).to eql(true)
          expect(::File.exists?(output_files[1])).to eql(true)

          ::File.write(output_files[0], "test data")
          ::File.write(output_files[1], "test data")

          expect(::File.read(output_files[0])).to_not eql("")
          expect(::File.read(output_files[1])).to_not eql("")

          raise "Oh no!"
        end
      end.to raise_error("Oh no!")

      built_paths = ['a', 'b'].map { |output_filename| storage.data_path(project_name: 'project', cell_hash: 'ch', data_filename: output_filename) }
      expect(::File.exists?(built_paths[0])).to eql(false)
      expect(::File.exists?(built_paths[1])).to eql(false)

      storage.with_build_transaction(project_name: 'project', ch: 'ch', output_filenames: ['a', 'b']) do |output_files|
        output_files = output_files.map { |of| of.data_path(storage) }
        expect(::File.exists?(output_files[0])).to eql(true)
        expect(::File.exists?(output_files[1])).to eql(true)

        expect(::File.read(output_files[0])).to eql("")
        expect(::File.read(output_files[1])).to eql("")

        ::File.write(output_files[0], "stuff 1")
        ::File.write(output_files[1], "stuff 2")
      end

      expect(::File.exists?(built_paths[0])).to eql(true)
      expect(::File.exists?(built_paths[1])).to eql(true)

      expect(::File.read(built_paths[0])).to eql("stuff 1")
      expect(::File.read(built_paths[1])).to eql("stuff 2")
    end
  end

  describe '.cell_hash' do
    before(:each) do
      @input_file_a = Vulcan::Storage::StorageFile.new(project_name: 'project', cell_hash: 'abc', data_filename: 'r', logical_name: 'a', prefix: 'a')
      @input_file_b = Vulcan::Storage::StorageFile.new(project_name: 'project', cell_hash: 'abc', data_filename: 'r', logical_name: 'b', prefix: 'a')
      @input_files = [@input_file_a, @input_file_b]
      @output_filenames = ['a', 'b']
      @session_key = 'session_key'
      @script = 'script'
      @raw_hash = nil
      @project_name = 'project'
    end

    def hash
      Vulcan::Storage.cell_hash(project_name: @project_name, input_files: @input_files, output_filenames: @output_filenames, session_key: @session_key, script: @script, raw_hash: @raw_hash)
    end

    it 'is stable for input file ordering' do
      expect { @input_files.reverse! }.to_not change { hash }
    end

    it 'is stable for input file prefix' do
      expect { @input_file_a.prefix += 'moo' }.to_not change { hash }
    end

    it 'is stable for output file name ordering' do
      expect { @output_filenames.reverse! }.to_not change { hash }
    end

    it 'changes by the input file hashes' do
      expect { @input_file_a.cell_hash += 'b' }.to change { hash }
    end

    it 'changes when the input is a script or a raw_hash' do
      expect { @raw_hash = @script; @script = nil }.to change { hash }
    end

    it 'changes when session_key changes' do
      expect { @session_key += 'session_key_more' }.to change { hash }
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
      expect(hash(Vulcan::Storage::StorageFile.new(**values))).to_not eql('')
      expect(hash(Vulcan::Storage::StorageFile.new(**values))).to eql(hash(values))
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
