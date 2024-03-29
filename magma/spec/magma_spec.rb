require_relative '../lib/magma'
require 'yaml'

def disconnect_attribute(model, att_name)
  model.attributes.delete(att_name)
end

def reconnect_attribute model, att_name, value
  model.attributes[att_name] = value
end

def reset_models project
  Magma.instance.get_project(project).instance_variable_set("@models",nil)
end

describe Magma do
  describe '#get_model' do
    it 'returns the class for a given model name' do
      model = Magma.instance.get_model('labors','project')

      expect(model).to eq(Labors::Project)
    end
    it 'raises if the model does not exist' do
      expect {
        Magma.instance.get_model('labors','nonexistent_model')
      }.to raise_error(NameError)
    end
  end

  describe '#validate_models' do
    context 'project model' do
      before(:each) do
        # Delete the constant
        @project = Labors::Project
        Labors.send(:remove_const,:Project)

        # Reset the Magma::Project object's model cache
        reset_models(:labors)

        # Unlink the labor model's project attribute
        @project_att = disconnect_attribute(Labors::Labor,:project)
      end
      after(:each) do
        Labors::Project = @project

        reset_models(:labors)

        reconnect_attribute(Labors::Labor, :project, @project_att)
      end
      it 'complains if no project model is present' do
        expect{Magma.instance.validate_models}.to raise_error(Magma::ValidationError)
      end
    end

    context 'orphan models' do
      before(:each) do
        # Unlink the labor model and prize model
        @labor_att = disconnect_attribute(Labors::Prize,:labor)
        @prize_att = disconnect_attribute(Labors::Labor,:prize)
      end
      after(:each) do
        reconnect_attribute(Labors::Labor, :prize, @prize_att)
        reconnect_attribute(Labors::Prize, :labor, @labor_att)
      end
      it 'complains if there are orphan models' do
        expect{Magma.instance.validate_models}.to raise_error(Magma::ValidationError)
      end
    end

    it "raises an error if a model has an attribute that isn't backed by a column in the model's table" do
      Labors::Monster.attributes[:size] = Magma::IntegerAttribute.new(
        attribute_name: "size",
        magma_model: Labors::Monster
      )

      expect{ Magma.instance.validate_models }.to raise_error(Magma::ValidationError)

      Labors::Monster.attributes.delete(:size)
    end
  end

  context 'hold file' do
    it 'writes the hold file' do
      Magma.instance.write_hold_file

      expect(::File.exists?(Magma.instance.config(:hold_file))).to be_truthy

      ::File.unlink(Magma.instance.config(:hold_file))
    end

    it 'removes the hold file' do
      Magma.instance.write_hold_file

      Magma.instance.remove_hold_file

      expect(::File.exists?(Magma.instance.config(:hold_file))).to be_falsy
    end

    it 'passes the healthcheck if hold file exists and is not expired' do
      Magma.instance.write_hold_file

      status = system("bin/healthcheck")

      ::File.unlink(Magma.instance.config(:hold_file))

      expect(status).to be_truthy
    end

    it 'runs the healthcheck if hold file does not exist' do
      status = system("bin/healthcheck")

      expect(status).to be_falsy
    end

    it 'runs the healthcheck if hold file is expired' do
      Timecop.freeze(DateTime.now-Magma.instance.config(:hold_interval)-1)
      Magma.instance.write_hold_file
      Timecop.return

      status = system("bin/healthcheck")

      ::File.unlink(Magma.instance.config(:hold_file))

      expect(status).to be_falsy
    end
  end
end
