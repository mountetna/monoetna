describe Polyphemus::RunEtlJob do
  let(:command) { Polyphemus.instance.subcommands['run_etl_job'] }

  context 'running etl jobs' do
    before(:all) do
      Timecop.freeze
      class Polyphemus
        class DummyJob < Polyphemus::Job
          def self.as_json
            {
              name: "dummy",
              schema: {
                type: "object",
                properties: {
                  foo: { type: "integer" },
                  bar: { enum: [ "baz", "qux" ] }
                },
                required: ["foo"],
                additionalProperties: false
              }
            }
          end

          def run
            puts "Here is some output"
          end
        end
      end
    end

    def create_dummy_etl(opts)
      create(:etl_config, {project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy"}.merge(opts))
    end

    after(:all) do
      Polyphemus::Job.list.delete(Polyphemus::DummyJob)
      Polyphemus.send(:remove_const,:DummyJob)
      Timecop.return
    end

    it 'runs a one-time Etl job' do
      etl = create_dummy_etl(etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_ONCE)
      
      command.execute

      etl.refresh

      expect(etl.ran_at).to be_within(1).of(Time.now)
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.status).to eq(Polyphemus::EtlConfig::STATUS_COMPLETED)
      expect(etl.output).to eq("Here is some output\n")
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
