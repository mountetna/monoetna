describe 'Magma Commands' do
  let(:magma_instance) { double('Magma') }

  describe Magma::Migrate do
    subject(:global_migrate) { described_class.new.execute(version: version) }
    let(:directory) { "db/migrations" }
    let(:table) {"schema_info_spec_labors" }

    before do
      Sequel.extension(:migration)
      allow(Magma.instance).to receive(:db).and_return(magma_instance)
    end

    describe 'when a version is specified' do
      let(:version) { '001' }

      before do
        allow(Sequel::Migrator).to receive(:run).with(magma_instance, directory, target: version.to_i).and_return(true)
      end

      it 'calls run with a version number' do
        global_migrate

        expect(Sequel::Migrator)
          .to have_received(:run)
          .with(magma_instance, directory, target: version.to_i)
          .once
      end
    end

    describe 'without a version specified' do
      let(:version) { nil }

      before do
        allow(Sequel::Migrator).to receive(:run).with(magma_instance, directory).and_return(true)
      end

      it 'calls run without a version number' do
        global_migrate

        expect(Sequel::Migrator)
          .to have_received(:run)
          .with(magma_instance, directory)
          .once
      end
    end
  end

  describe Magma::Console do
    subject(:console) { described_class.new.execute }
    before do
      require 'irb'
      allow(ARGV).to receive(:clear)
      allow(IRB).to receive(:start)
    end

    it 'calls ARGV and IRB' do
      console

      expect(ARGV).to have_received(:clear).once
      expect(IRB).to have_received(:start).once
    end
  end

  describe Magma::Unload do
    let(:project_name) { 'project_name' }
    let(:model_name) { 'model_name' }
    let(:tsv_writer_instance) { double('tsv_writer') }
    let(:model_double) { double('model') }
    let(:retrieval_double) { double('retrieval') }
    let(:payload_double) { double('payload') }

    subject(:unload) { described_class.new.execute(project_name, model_name) }

    before do
      allow(Magma).to receive(:instance).and_return(magma_instance)
      allow(magma_instance).to receive(:get_model).with(project_name, model_name).and_return(model_double)
      allow(Magma::Retrieval).to receive(:new).with(model_double, 'all', 'all', page: 1, page_size: 100_000).and_return(retrieval_double)
      allow(Magma::Payload).to receive(:new).and_return(payload_double)
      allow(Magma::TSVWriter).to receive(:new).with(model_double, retrieval_double, payload_double).and_return(tsv_writer_instance)
      allow(tsv_writer_instance).to receive(:write_tsv)
    end

    it 'calls TSVWriter.new once' do
      unload

      expect(tsv_writer_instance).to have_received(:write_tsv).once
    end
  end

 describe Magma::Load do
    subject(:loader_execute) { described_class.new.execute('loader name', 'args') }

    let(:loader) { double('loader', loader_name: 'loader name', description: 'description') }
    let(:loaders_array) { [loader] }
    let(:expected) { 'yup' }

    before do
      allow(Magma).to receive(:instance).and_return(magma_instance)
      allow(magma_instance).to receive(:find_descendents).with(Magma::Loader).and_return(loaders_array)
      allow(loader).to receive(:new).and_return(loader)
      allow(loader).to receive(:load).with('args')
      allow(loader).to receive(:dispatch)
    end

    it 'displays available loaders' do
      loader_execute

      expect(loader).to have_received(:dispatch).once
    end
  end

  describe Magma::AutoDisconnect do
    let(:project_name) { 'labors' }
    let(:project) { create(:project, name: project_name) }

    before(:each) do
      @labor1 = create(:labor, name: 'Nemean Lion', project: project)
      @labor2 = create(:labor, name: 'Lernean Hydra', project: project)
      
      @monster1 = create(:monster, name: 'Lion', labor: @labor1)
      @monster2 = create(:monster, name: 'Hydra', labor: @labor2)
      @monster3 = create(:monster, name: 'Disconnected Monster') # Already disconnected
      
      @victim1 = create(:victim, name: 'Villager1', monster: @monster1)
      @victim2 = create(:victim, name: 'Villager2', monster: @monster2)
    end

    describe '#execute without --auto flag' do
      it 'disconnects all records of a single model' do
        expect(Labors::Monster.where(Sequel.~(labor_id: nil)).count).to eq(2)
        expect(Labors::Monster.where(labor_id: nil).count).to eq(1)

        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        expect {
          command.execute(project_name, 'monster', auto: false)
        }.to output(/Successfully disconnected 2 record\(s\)/).to_stdout

        expect(Labors::Monster.where(labor_id: nil).count).to eq(3)
        
        @monster1.refresh
        @monster2.refresh
        expect(@monster1.labor).to be_nil
        expect(@monster2.labor).to be_nil
      end

      it 'reports when no connected records are found' do
        # First disconnect all monsters
        Labors::Monster.where(Sequel.~(labor_id: nil)).update(labor_id: nil)

        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        expect {
          command.execute(project_name, 'monster', auto: false)
        }.to output(/No connected records found/).to_stdout
      end

      it 'reports error for models without parent attribute' do
        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        expect {
          command.execute(project_name, 'project', auto: false)
        }.to output(/does not have a parent attribute/).to_stdout
      end

      it 'preserves child relationships when disconnecting parent' do
        # Disconnect monsters
        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        command.execute(project_name, 'monster', auto: false)

        # Victims should still be connected to their monsters
        @victim1.refresh
        @victim2.refresh
        expect(@victim1.monster).to eq(@monster1)
        expect(@victim2.monster).to eq(@monster2)
      end
    end

    describe '#execute with --auto flag' do
      it 'disconnects records up the hierarchy' do
        # Create a deeper hierarchy: labor -> monster -> victim
        expect(Labors::Victim.where(Sequel.~(monster_id: nil)).count).to eq(2)
        expect(Labors::Monster.where(Sequel.~(labor_id: nil)).count).to eq(2)

        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        expect {
          command.execute(project_name, 'victim', auto: true)
        }.to output(/Auto mode.*victim.*monster.*labor/m).to_stdout

        # All victims should be disconnected
        expect(Labors::Victim.where(monster_id: nil).count).to eq(2)
        
        # All monsters should be disconnected
        expect(Labors::Monster.where(labor_id: nil).count).to eq(3)

        @victim1.refresh
        @victim2.refresh
        @monster1.refresh
        @monster2.refresh

        expect(@victim1.monster).to be_nil
        expect(@victim2.monster).to be_nil
        expect(@monster1.labor).to be_nil
        expect(@monster2.labor).to be_nil
      end

      it 'stops at project level and does not disconnect project' do
        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        output = nil
        expect {
          command.execute(project_name, 'monster', auto: true)
        }.to output { |out| output = out }.to_stdout

        # Should not mention project in the hierarchy
        expect(output).not_to match(/\bproject\b/)
        expect(output).to match(/monster/)
        expect(output).to match(/labor/)
      end

      it 'disconnects in correct order (children before parents)' do
        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        disconnection_order = []
        
        # Track when disconnect_model is called
        original_disconnect = command.method(:disconnect_model)
        allow(command).to receive(:disconnect_model) do |model|
          disconnection_order << model.model_name.to_s
          original_disconnect.call(model)
        end

        command.execute(project_name, 'victim', auto: true)

        # Victims should be disconnected before monsters, monsters before labors
        victim_index = disconnection_order.index('victim')
        monster_index = disconnection_order.index('monster')
        labor_index = disconnection_order.index('labor')

        expect(victim_index).to be < monster_index if victim_index && monster_index
        expect(monster_index).to be < labor_index if monster_index && labor_index
      end
    end

    describe 'edge cases' do
      it 'handles models with no records gracefully' do
        # Create a new labor with no monsters
        labor3 = create(:labor, name: 'Empty Labor', project: project)

        # Delete all monsters first
        Labors::Monster.where(Sequel.~(labor_id: nil)).delete

        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        
        expect {
          command.execute(project_name, 'monster', auto: false)
        }.to output(/No connected records found/).to_stdout
      end

      it 'handles table-type models' do
        # Prize is a table-type model
        prize1 = create(:prize, labor: @labor1, name: 'Gold', worth: 100)
        prize2 = create(:prize, labor: @labor2, name: 'Silver', worth: 50)

        expect(Labors::Prize.where(Sequel.~(labor_id: nil)).count).to eq(2)

        command = described_class.new
        command.setup(YAML.load(File.read('config.yml')))
        command.execute(project_name, 'prize', auto: false)

        expect(Labors::Prize.where(labor_id: nil).count).to eq(2)
        
        prize1.refresh
        prize2.refresh
        expect(prize1.labor_id).to be_nil
        expect(prize2.labor_id).to be_nil
      end
    end
  end
end
