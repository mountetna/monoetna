require_relative '../lib/commands'

describe Etna::GenerateCompletionScript do
  def test_completion(*parts)
    parts = ['etna'] + parts.map { |p| p == '' ? '" "' : p }
    cmd, _, _ = EtnaApp.instance.find_command('generate_completion_script')
    cmd.define_singleton_method(:program_name) { 'etna' }
    cmd.execute
    `bash -c 'source #{__dir__}/../etna.completion && #{__dir__}/get_completions.sh #{parts.join(' ')}'`.split("\n")
  end

  subject { test_completion(*args) }

  describe 'root level flag completions' do
    let(:args) { [''] }
    it { should include('--environment') }
    it { should include('config') }
    it { should include('help') }

    describe 'after another command' do
      let(:args) { ['config', ''] }
      it { should include('--environment') }
      it { should include('help') }
      it { should include('set') }

      describe 'completed with argument' do
        let(:args) { [ 'config', '--environment', 'development', ''] }
        it { should_not include('--environment') }
        it { should include('set') }
        it { should include('help') }
      end
    end

    describe 'with string flag completion' do
      let(:args) { ['--environment', ''] }
      it { should include('development') }
      it { should_not include('config') }
      it { should_not include('--environment') }

      describe 'completed with argument' do
        let(:args) { ['--environment', 'development', ''] }
        it { should_not include('--environment') }
        it { should include('config') }
        it { should include('help') }
      end
    end
  end

  describe 'command level completions' do
    let(:args) { ['config', 'set', ''] }
    it { should include('__host__') }

    describe 'after argument completion' do
      let(:args) { ['config', 'set', '__host__', ''] }
      it { should include('--environment') }
      it { should include('--ignore-ssl') }
      it { should_not include('__host__') }
      it { should_not include('help') }

      describe 'with a boolean flag completion' do
        let(:args) { ['config', 'set', '__host__', '--ignore-ssl', ''] }
        it { should include('--environment') }
        it { should_not include('--ignore-ssl') }
        it { should_not include('__host__') }
        it { should_not include('help') }

        describe 'and completed all remaining flags' do
          let(:args) { ['config', 'set', '__host__', '--ignore-ssl', '--environment', 'development', '' ] }
          it { should eql(['DEFAULT']) }
        end
      end
    end
  end
end