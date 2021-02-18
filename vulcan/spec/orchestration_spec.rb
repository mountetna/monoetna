describe Vulcan::Orchestration do
  let(:storage) { Vulcan::Storage.new }
  let(:inputs) { {} }
  let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs) }
  let(:orchestration) { session.orchestration }

  before(:each) do
    FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)
  end

  def primary_outputs
    orchestration.build_target_for(:primary_outputs)
  end

  def unique_paths
    orchestration.unique_paths
  end

  def build_targets
    orchestration.build_targets_for_paths
  end

  def should_builds
    build_targets.map { |b| b.map { |bt| bt.should_build?(storage) } }
  end

  def cell_hashes
    build_targets.map { |b| b.map { |bt| bt.cell_hash } }
  end

  def next_buildable
    orchestration.next_runnable_build_targets(storage).first
  end

  def cell_hashes_same(old_cell_hashes)
    cell_hashes.zip(old_cell_hashes).map do |a, b|
      a.zip(b).map do |at, bt|
        at == bt
      end
    end
  end

  describe '#unique_paths' do
    it 'works' do
      expect(unique_paths).to eql([
          [:primary_inputs, "firstAdd", "finalStep", :primary_outputs],
          [:primary_inputs, "firstAdd", "pickANum", "finalStep"]
      ])
    end
  end

  describe '#build_targets_for_paths' do
    describe 'should_builds / next_runnable_build_targets' do
      it 'works' do
        orchestration.load_json_inputs!(storage)

        expect(should_builds).to eql([
            [false, false, false, false],
            [false, false, false, false],
        ])
        expect(orchestration.next_runnable_build_targets(storage)).to eql([])

        # supplying at least the non default value should move things along
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        orchestration.load_json_inputs!(storage)

        expect(should_builds).to eql([
            [true, false, false, false],
            [true, false, false, false],
        ])

        expect(orchestration.next_runnable_build_targets(storage).map(&:cell_hash)).to eql([
            orchestration.build_target_for(:primary_inputs).cell_hash
        ])

        expect(next_buildable.script).to be_a(Hash)
        orchestration.run!(storage: storage, build_target: next_buildable)

        expect(should_builds).to eql([
            [false, true, false, false],
            [false, true, false, false],
        ])
      end
    end

    describe 'cell_hash calculations' do
      it 'works' do
        expect(cell_hashes.flatten.sort.uniq.length).to eql(5)

        prev_cell_hashes = cell_hashes
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [false, false, false, false],
            [false, false, false, false]
        ])

        # Test verification: check the falsifiable condition that the reason these tests pass is because this function
        # does not work.
        prev_cell_hashes = cell_hashes
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true]
        ])

        prev_cell_hashes = cell_hashes
        session.define_user_input(["primary_inputs", "notARealInput"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true]
        ])

        prev_cell_hashes = cell_hashes
        session.define_user_input([:primary_inputs, "alsoNotAnInput"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true]
        ])

        # Defining an input only effects downstream cells
        prev_cell_hashes = cell_hashes
        session.define_user_input(["pickANum", "num"], 543)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, false, false],
            [true, true, false, false]
        ])
      end
    end
  end

  describe 'e2e' do
    it 'works' do
      expect(orchestration.run_until_done!(storage).length).to eql(0)
      expect(primary_outputs.is_built?(storage)).to eql(false)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      expect(orchestration.run_until_done!(storage).length).to eql(2)
      expect(primary_outputs.is_built?(storage)).to eql(false)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      session.define_user_input(["pickANum", "num"], 543)
      expect(orchestration.run_until_done!(storage).length).to eql(3)
      expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
    end
  end
end
