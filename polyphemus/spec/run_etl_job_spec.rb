describe Polyphemus::RunEtlJob do
  let(:command) { Polyphemus.instance.subcommands['run_etl_job'] }

  context 'running etl jobs' do
    before(:all) do
      create_dummy_job
    end

    before(:each) do
      Timecop.freeze
      stub_task_token
    end

    after(:each) do
      Timecop.return
    end

    after(:all) do
      remove_dummy_job
    end

    it 'runs a one-time Etl job' do
      etl = create_dummy_etl(etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_ONCE, params: { problem: 'absent', whippit: true } )
      
      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.status).to eq(Polyphemus::EtlConfig::STATUS_COMPLETED)
      expect(etl.params).to eq({})
      expect(etl.output).to eq("Here is some output\n")
    end

    it 'catches errors' do
      etl = create_dummy_etl(etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_ONCE, params: { problem: 'absent', whippit: true } )

      Polyphemus::DummyJob.define_method :run do
        raise 'Uh-oh'
      end
      
      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.status).to eq(Polyphemus::EtlConfig::STATUS_ERROR)
      expect(etl.params).to eq('problem' => 'absent', 'whippit' => true)
      expect(etl.output).to eq("\nUh-oh")

      Polyphemus::DummyJob.define_method :run do
        puts 'Here is some output'
      end
    end

    it 'runs an interval Etl job' do
      etl = create_dummy_etl(etl: "dummy", run_interval: 60)
      
      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.run_interval).to eq(60)
      expect(etl.status).to eq(Polyphemus::EtlConfig::STATUS_COMPLETED)
      expect(etl.output).to eq("Here is some output\n")
    end

    it 'waits to run an interval Etl job again' do
      etl = create_dummy_etl(etl: "dummy", run_interval: 60, ran_at: DateTime.now, output: "It did not run")
      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.output).to eq("It did not run")
    end

    it 'runs an interval Etl job again after waiting' do
      etl = create_dummy_etl(etl: "dummy", run_interval: 60, ran_at: DateTime.now - 80, output: "It did not run")

      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.output).to eq("Here is some output\n")
    end

    it 'ignores an Etl job which should never run' do
      etl = create_dummy_etl(etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      
      command.execute

      etl.refresh

      expect(etl.ran_at).to eq(nil)
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.status).to eq(nil)
      expect(etl.output).to eq(nil)
    end

    it 'skips an archived job' do
      etl = create_dummy_etl(etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_ONCE, archived: true)

      command.execute

      etl.refresh

      expect(etl.ran_at).to eq(nil)
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_ONCE)
    end
  end
end
