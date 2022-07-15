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

  def serialized_step_path
    orchestration.serialized_step_path
  end

  def build_targets
    orchestration.build_targets_for_paths
  end

  def should_builds
    build_targets.map { |b| b.map { |step_name, bt| bt.should_build?(storage) } }
  end

  def cell_hashes
    build_targets.map { |b| b.map { |step_name,bt| bt.cell_hash } }
  end

  def next_buildable
    orchestration.next_runnable_build_targets(storage).first.last
  end

  def cell_hashes_same(old_cell_hashes)
    cell_hashes.zip(old_cell_hashes).map do |a, b|
      a.zip(b).map do |at, bt|
        at == bt
      end
    end
  end

  describe '#serialized_step_path' do
    it 'works' do
      expect(serialized_step_path).to eql([
          [:primary_inputs_someIntWithoutDefault, :primary_inputs_someInt, "firstAdd", "pickANum", "finalStep", "aPlot", :primary_outputs],
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
            [false, false, false, false],
            [false, false]
        ])
        expect(orchestration.next_runnable_build_targets(storage)).to eql([])

        # supplying at least the non default value should move things along
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        orchestration.load_json_inputs!(storage)

        expect(should_builds).to eql([
            [false, true, false, false],
            [false, true, false, false],
            [false, true, false, false],
            [false, true]
        ])

        expect(orchestration.next_runnable_build_targets(storage).map(&:last).map(&:cell_hash)).to eql([
            orchestration.build_target_for('firstAdd').cell_hash
        ])

        orchestration.run!(
          storage: storage,
          build_target: next_buildable,
          token: 'test-token')

        session.define_user_input(['pickANum', "num"], 456)
        orchestration.load_json_inputs!(storage)

        expect(should_builds).to eql([
            [false, false, false, false],
            [false, false, false, false],
            [false, false, true, false],
            [false, false]
        ])
      end
    end

    describe 'cell_hash calculations' do
      it 'works' do
        expect(cell_hashes.flatten.sort.uniq.length).to eql(7)


        prev_cell_hashes = cell_hashes
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [false, false, false, false],
            [false, false, false, false],
            [false, false, false, false],
            [true, false],
        ])

        # Test verification: check the falsifiable condition that the reason these tests pass is because this function
        # does not work.
        prev_cell_hashes = cell_hashes
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true],
            [true, true, true, true],
            [true, true],
        ])

        prev_cell_hashes = cell_hashes
        session.define_user_input(["primary_inputs", "notARealInput"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true],
            [true, true, true, true],
            [true, true]
        ])

        prev_cell_hashes = cell_hashes
        session.define_user_input([:primary_inputs, "alsoNotAnInput"], 123)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, true, true],
            [true, true, true, true],
            [true, true, true, true],
            [true, true]
        ])

        # Defining an input only effects downstream cells
        prev_cell_hashes = cell_hashes
        session.define_user_input(["pickANum", "num"], 543)
        expect(cell_hashes_same(prev_cell_hashes)).to eql([
            [true, true, false, false],
            [true, true, false, false],
            [true, true, false, false],
            [true, true]
        ])
      end
    end
  end

  describe 'e2e' do
    it 'works' do
      expect(orchestration.run_until_done!(storage).length).to eql(0)
      expect(primary_outputs.is_built?(storage)).to eql(false)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      expect(orchestration.run_until_done!(storage).length).to eql(1)
      expect(primary_outputs.is_built?(storage)).to eql(false)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      session.define_user_input(["pickANum", "num"], 543)
      expect(orchestration.run_until_done!(storage).length).to eql(3)
      expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
    end

    it 'reports cell errors correctly' do
      expect(orchestration.run_until_done!(storage).length).to eql(0)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 'abc-not-an-int')

      expect {
        orchestration.run_until_done!(storage)
      }.to raise_error(Vulcan::Orchestration::RunErrors)
    end

    describe 'specifying an archimedes sha' do
      let(:figure) {
        create_figure(
          workflow_name: "test_workflow.cwl",
          dependencies: {
            "etnaagent/archimedes": "sha256:0336c5101c0a089a2bb38d2fa01c0747f4f6bd615e0326a567a256e8aa04d4b0"
          }
        )
      }
      let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
    
      it 'works' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(orchestration.run_until_done!(storage).length).to eql(1)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        session.define_user_input(["pickANum", "num"], 543)
        expect(orchestration.run_until_done!(storage).length).to eql(3)
        expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
      end

      it 'reports cell errors correctly' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 'abc-not-an-int')

        expect {
          orchestration.run_until_done!(storage)
        }.to raise_error(Vulcan::Orchestration::RunErrors)
      end

      describe 'when dependencies are ignored' do
        it 'generates new cell hashes' do
          expect(orchestration.run_until_done!(storage).length).to eql(0)
          expect(primary_outputs.is_built?(storage)).to eql(false)
          session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          expect(orchestration.run_until_done!(storage).length).to eql(1)
          expect(primary_outputs.is_built?(storage)).to eql(false)
          session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          session.define_user_input(["pickANum", "num"], 543)
          expect(orchestration.run_until_done!(storage).length).to eql(3)

          final_data_path_for_session = primary_outputs.build_outputs['the_result'].data_path(storage)
          
          ignored_session = Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', {}, reference_figure_id: figure.id, ignore_dependencies: true)
          ignored_orchestration = ignored_session.orchestration

          expect(ignored_orchestration.run_until_done!(storage).length).to eql(0)

          ignored_primary_outputs = ignored_orchestration.build_target_for(:primary_outputs)
          expect(ignored_primary_outputs.is_built?(storage)).to eql(false)
          ignored_session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          expect(ignored_orchestration.run_until_done!(storage).length).to eql(1)
          expect(ignored_primary_outputs.is_built?(storage)).to eql(false)
          ignored_session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          ignored_session.define_user_input(["pickANum", "num"], 543)
          expect(ignored_orchestration.run_until_done!(storage).length).to eql(3)

          final_data_path_for_ignored_session = ignored_primary_outputs.build_outputs['the_result'].data_path(storage)

          # If the cell_hash / result path has changed, then we'll assume
          #   the ignore_dependencies flag worked ... at the very least,
          #   the cache responded correctly. 
          expect(final_data_path_for_ignored_session).not_to eq(final_data_path_for_session)
        end
      end
    end

    describe 'specifying a non-existent archimedes sha' do
      let(:figure) {
        create_figure(
          workflow_name: "test_workflow.cwl",
          dependencies: {
            "etnaagent/archimedes": "sha256:f3b3b28a45160805bb16542c9531888519430e9e6d6ffc09d72261b0d26ff74f"
          }
        )
      }
      let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
    
      it 'works' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(orchestration.run_until_done!(storage).length).to eql(1)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        session.define_user_input(["pickANum", "num"], 543)
        expect(orchestration.run_until_done!(storage).length).to eql(3)
        expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
      end

      it 'reports cell errors correctly' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 'abc-not-an-int')

        expect {
          orchestration.run_until_done!(storage)
        }.to raise_error(Vulcan::Orchestration::RunErrors)
      end
    end

    describe 'with an empty archimedes sha' do
      let(:figure) {
        create_figure(
          workflow_name: "test_workflow.cwl",
          dependencies: {
            "etnaagent/archimedes": ""
          }
        )
      }
      let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
    
      it 'works' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(orchestration.run_until_done!(storage).length).to eql(1)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        session.define_user_input(["pickANum", "num"], 543)
        expect(orchestration.run_until_done!(storage).length).to eql(3)
        expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
      end

      it 'reports cell errors correctly' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 'abc-not-an-int')

        expect {
          orchestration.run_until_done!(storage)
        }.to raise_error(Vulcan::Orchestration::RunErrors)
      end
    end

    describe 'node' do
      let(:figure) {
        create_figure(
          workflow_name: "test_node_workflow.cwl",
          dependencies: {
            "etnaagent/archimedes-node": "sha256:969d1df72b9a6f78a17600f78316d0a48b94bb18679f9df04c7d28629891536a",
            "etnaagent/archimedes": "sha256:0336c5101c0a089a2bb38d2fa01c0747f4f6bd615e0326a567a256e8aa04d4b0"
          }
        )
      }
      let(:session) { Session.new_session_for('project', 'test_node_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
    
      it 'works' do
        expect(orchestration.run_until_done!(storage).length).to eql(2)
        expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql({
          'col1' => {
            0 => 1,
            1 => 2
          },
          'col2' => {
            0 => 'abc',
            1 => 'xyz'
          }
        }.to_json)
      end
    end

    describe 'using snapshot scripts instead of current scripts' do
      let(:figure) {
        create_figure(
          workflow_name: "test_workflow.cwl"
        )
      }
      let(:session) { Session.new_session_for('project', 'test_workflow.cwl', 'storage_key', inputs, reference_figure_id: figure.id) }
    
      before(:each) do
        # Inject an old "add" script that actually subtracts
        old_scripts = figure.workflow_snapshot.scripts
        figure.workflow_snapshot.update(
          scripts: {
            "aPlot"=>nil,
            "finalStep"=>
             "from archimedes.functions.dataflow import output_path, input_path\n\na = int(open(input_path('a'), 'r').read())\nb = int(open(input_path('b'), 'r').read())\n\nwith open(output_path('sum'), 'w') as output_file:\n    output_file.write(str(a - b))",
            "firstAdd"=>
             "from archimedes.functions.dataflow import output_path, input_path\n\na = int(open(input_path('a'), 'r').read())\nb = int(open(input_path('b'), 'r').read())\n\nwith open(output_path('sum'), 'w') as output_file:\n    output_file.write(str(a - b))",
            "pickANum"=>nil
          }
        )
      end

      it 'works' do
        expect(orchestration.run_until_done!(storage).length).to eql(0)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        expect(orchestration.run_until_done!(storage).length).to eql(1)
        expect(primary_outputs.is_built?(storage)).to eql(false)
        session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
        session.define_user_input(["pickANum", "num"], 543)
        expect(orchestration.run_until_done!(storage).length).to eql(3)
        expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("-466") # 200 - 123 - 543
      end

      describe 'falls back to current script if not stored in snapshot' do
        before(:each) do
          old_scripts = figure.workflow_snapshot.scripts
          figure.workflow_snapshot.update(
            scripts: {
              "aPlot"=>nil,
              "firstAdd"=>
               "from archimedes.functions.dataflow import output_path, input_path\n\na = int(open(input_path('a'), 'r').read())\nb = int(open(input_path('b'), 'r').read())\n\nwith open(output_path('sum'), 'w') as output_file:\n    output_file.write(str(a - b))",
              "pickANum"=>nil
            }
          )
        end
  
        it 'works' do
          expect(orchestration.run_until_done!(storage).length).to eql(0)
          expect(primary_outputs.is_built?(storage)).to eql(false)
          session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          expect(orchestration.run_until_done!(storage).length).to eql(1)
          expect(primary_outputs.is_built?(storage)).to eql(false)
          session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
          session.define_user_input(["pickANum", "num"], 543)
          expect(orchestration.run_until_done!(storage).length).to eql(3)
          expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("620") # 200 - 123 + 543
        end
      end
    end
  end
end
