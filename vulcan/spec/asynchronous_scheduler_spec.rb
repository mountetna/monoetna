describe Vulcan::AsynchronousScheduler do
  let(:storage) { Vulcan::Storage.new }
  let(:inputs) { {} }
  let(:session) { Session.new_session_for('project', 'test_concurrent_workflow.cwl', 'async-scheduler-test', inputs) }
  let(:orchestration) { session.orchestration }
  let(:scheduler) { orchestration.scheduler }

  before(:each) do
    FileUtils.rm_rf(storage.data_root) if ::File.exist?(storage.data_root)
  end

  def primary_outputs
    orchestration.build_target_for(:primary_outputs)
  end

  describe 'e2e' do
    it 'works' do
      expect(scheduler.schedule_more!(orchestration: orchestration, token: '', storage: storage)).to eql([])
      session.define_user_input([:primary_inputs, "a"], 123)
      session.define_user_input([:primary_inputs, "b"], 234)
      session.define_user_input(["pickANum", "num"], 51)
      session.define_user_input(["otherPickANumber", "num"], 51)

      orchestration.load_json_inputs!(storage)

      scheduler.instance_variable_get(:@_test_gating_semaphore).synchronize do
        expect(scheduler.schedule_more!(orchestration: orchestration, token: '', storage: storage).length).to eql(2)
        expect(scheduler.schedule_more!(orchestration: orchestration, token: '', storage: storage).length).to eql(0)
      end

      scheduler.join_all
      expect(::File.read(primary_outputs.build_outputs['the_result'].data_path(storage))).to eql("306")
    end

    it 'marks and reports status correctly' do
      expect(orchestration.run_until_done!(storage).length).to eql(0)
      session.define_user_input([:primary_inputs, "a"], 'not an int')
      session.define_user_input([:primary_inputs, "b"], 234)

      bt = orchestration.build_target_for('firstAdd')

      expect(scheduler.status(storage: storage, build_target: bt)).to eql({ status: 'pending' })

      expect {
        orchestration.load_json_inputs!(storage)
        scheduler.schedule_more!(orchestration: orchestration, token: '', storage: storage)
        scheduler.join_all
      }.to change {
        scheduler.status(storage: storage, build_target: bt)[:status]
      }.to('error')

      expect(scheduler.status(storage: storage, build_target: bt)[:error]).to match(
           /Python error while executing script: Traceback \(most recent call last\)/
      )

      session.define_user_input([:primary_inputs, "a"], 483)
      bt = orchestration.build_target_for('firstAdd')

      expect {
        orchestration.load_json_inputs!(storage)
        scheduler.schedule_more!(orchestration: orchestration, token: '', storage: storage)
        scheduler.join_all
      }.to change {
        scheduler.status(storage: storage, build_target: bt)[:status]
      }.to('complete')
    end
  end
end
