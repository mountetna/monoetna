describe Vulcan::AsynchronousScheduler do
  let(:storage) { Vulcan::Storage.new }
  let(:inputs) { {} }
  let(:session) { Session.new_session_for('project', 'test_concurrent_workflow.cwl', 'async-scheduler-test', inputs) }
  let(:orchestration) { session.orchestration }
  let(:scheduler) { orchestration.scheduler }

  def primary_outputs
    orchestration.build_target_for(:primary_outputs)
  end

  describe 'e2e' do
    it 'works' do
      expect(scheduler.schedule_more!(token: '', storage: storage)).to eql(0)
      session.define_user_input([:primary_inputs, "a"], 123)
      session.define_user_input([:primary_inputs, "b"], 234)
      session.define_user_input(["pickANum", "num"], 51)
      session.define_user_input(["otherPickANumber", "num"], 51)

      orchestration.load_json_inputs!(storage)
      p(orchestration.unique_paths)
      p(orchestration.next_runnable_build_targets(storage))
      expect(scheduler.schedule_more!(token: '', storage: storage)).to eql(2)
      expect(scheduler.schedule_more!(token: '', storage: storage)).to eql(0)
      scheduler.join_all

      # expect(orchestration.run_until_done!(storage).length).to eql(0)
      # expect(primary_outputs.is_built?(storage)).to eql(false)
      # session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      # expect(orchestration.run_until_done!(storage).length).to eql(2)
      # expect(primary_outputs.is_built?(storage)).to eql(false)
      # session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 123)
      # session.define_user_input(["pickANum", "num"], 543)
      # expect(orchestration.run_until_done!(storage).length).to eql(3)
      # expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("866")
    end

    it 'reports cell errors correctly' do
      expect(orchestration.run_until_done!(storage).length).to eql(0)
      session.define_user_input([:primary_inputs, "someIntWithoutDefault"], 'abc-not-an-int')

      expect {
        orchestration.run_until_done!(storage)
      }.to raise_error(Vulcan::Orchestration::RunErrors)
    end
  end
end
