require 'singleton'

describe Etna::Command do
  class ExecutorSuper
    include Singleton
  end

  class Executor < ExecutorSuper
    include Etna::CommandExecutor

    boolean_flags << '--switch-1'
    boolean_flags << '--switch-2'

    class Thing < Etna::Command
      def execute(a, b, c, **kwds)
        puts "#{a} #{b} #{c}"
      end
    end

    class NotAThing
    end

    class Other
      include Etna::CommandExecutor

      string_flags << '--string-flag-1'

      class Inner < Etna::Command
        string_flags << '--string-flag-2'

        def execute(*args)
          puts "in here #{args}"
        end
      end
    end
  end

  def run_command(*args)
    cmd, cmd_args, cmd_kwds = Executor.instance.find_command(*args)
    cmd.execute(*cmd_args, **cmd_kwds)
  end

  let(:all_commands) do

  end

  before(:each) do
    @output = []
    ([Executor.instance] + Executor.instance.all_subcommands).each do | executor|
      allow(executor).to receive(:puts).and_wrap_original { |m, *args| @output.push(*args) }
    end
  end

  it 'should support top level and nested help commands automatically' do
    run_command("help")
    expect(@output.first).to include("usage:")
    expect(@output.first).to include("/rspec")

    expect(@output[3]).to include("thing")
    expect(@output[3]).to include("<a> <b> <c> [flags...]")

    expect(@output[1]).to include("help")

    expect(@output[2]).to include("other")

    @output = []
    run_command("other", "help")

    expect(@output[0]).to include("usage:")
    expect(@output[0]).to include("/rspec other")

    expect(@output[1]).to include("help")
    expect(@output[2]).to include("inner")
  end

  it 'should default to help on missing command' do
    run_command("other")
    expect(@output[0]).to include("usage:")
    expect(@output[0]).to include("/rspec other")

    expect(@output[1]).to include("help")
    expect(@output[2]).to include("inner")
  end

  it 'supports running commands, passing along arguments correctly' do
    run_command('thing', 'valuea', 'valueb', 'valuec')
    expect(@output[0]).to eql("valuea valueb valuec")

    @output = []

    run_command("other", "inner", "vardic", "args")
    expect(@output[0]).to eql('in here ["vardic", "args"]')
  end

  describe 'with flags' do
    it 'supports flags at executor level' do
      run_command("--switch-1", "other", "inner", "vardic", "args")
      expect(@output[0]).to eql('in here ["vardic", "args", {:switch_1=>true}]')
    end

    it 'supports flags positioned further down args' do
      run_command("other", "inner", "vardic", "args", "--switch-1")
      expect(@output[0]).to eql('in here ["vardic", "args", {:switch_1=>true}]')
    end

    it 'does not support flags positioned above own scope' do
      expect do
        run_command("--string-flag-2", "value", "other", "inner", "vardic", "args")
      end.to raise_error { |err| err.message =~ /rspec does not recognize flag --string-flag-2/ }
    end

    it 'supports nested contexts having their own flags' do
      run_command("--switch-2", "other", "inner", "vardic", "args", "--string-flag-2", "value")
      expect(@output[0]).to eql('in here ["vardic", "args", {:switch_2=>true, :string_flag_2=>"value"}]')
    end
  end
end