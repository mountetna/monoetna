require_relative '../lib/command_builder'

RSpec.describe Vulcan::CommandBuilder do
  let(:builder) { Vulcan::CommandBuilder.new }

  describe '#add' do
    it 'adds a command with arguments' do
      builder.add('echo', 'hello')
      expect(builder.to_s).to eq('echo hello')
    end
  end

  describe '#pipe_to' do
    it 'adds a pipe to another command' do
      builder.add('echo', 'hello').pipe_to('grep', 'h')
      expect(builder.to_s).to eq('echo hello | grep h')
    end
  end

  describe '#redirect_to' do
    it 'redirects output to a file' do
      builder.add('echo', 'hello').redirect_to('file.txt')
      expect(builder.to_s).to eq('echo hello > file.txt')
    end
  end

  describe '#append_to' do
    it 'appends output to a file' do
      builder.add('echo', 'hello').append_to('file.txt')
      expect(builder.to_s).to eq('echo hello >> file.txt')
    end
  end

  describe '#background' do
    it 'runs the command in the background' do
      builder.add('sleep', '10').background
      expect(builder.to_s).to eq('sleep 10 &')
    end
  end

  describe '#add_raw' do
    it 'adds a raw command with complex structure' do
      builder.add_raw('echo hello && echo world | grep w > output.txt')
      expect(builder.to_s).to eq('echo hello && echo world | grep w > output.txt')
    end
  end

  describe '#to_s' do
    it 'builds a complex command string' do
      builder.add('echo', 'hello').pipe_to('grep', 'h').redirect_to('file.txt').background
      expect(builder.to_s).to eq('echo hello | grep h > file.txt &')
    end
  end
end
