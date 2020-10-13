require 'singleton'

describe Etna::Command do
  class ExecutorSuper
    include Singleton
  end

  class Executor < ExecutorSuper
    include Etna::CommandExecutor

    class Thing < Etna::Command
      def execute(a, b, c)
        puts "#{a} #{b} #{c}"
      end
    end

    class NotAThing
    end

    class Other
      include Etna::CommandExecutor

      class Inner < Etna::Command
        def execute(*args)
          puts "in here #{args}"
        end
      end
    end
  end

  def run_command(*args)
    cmd, cmd_args = Executor.instance.find_command(*args)
    cmd.execute(*cmd_args)
  end

  before(:each) do
    @output = []
    [Executor.instance, Executor.instance.subcommands['thing'], Executor.instance.subcommands['other'], Executor.instance.subcommands['other'].subcommands['inner']].each do | executor|
      allow(executor).to receive(:puts).and_wrap_original { |m, *args| @output.push(*args) }
    end
  end

  it 'should support top level and nested help commands automatically' do
    run_command("help")
    expect(@output.first).to include("usage:")
    expect(@output.first).to include("/rspec")

    expect(@output[3]).to include("thing")
    expect(@output[3]).to include("<a> <b> <c>")

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
end