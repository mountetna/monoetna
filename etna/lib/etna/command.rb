require 'rollbar'

# Commands resolution works by starting from a class that includes CommandExecutor, and searches for
# INNER (belong to their scope) classes that include :find_command (generally subclasses of Etna::Command or
# classes that include CommandExecutor).
# A command completes the find_command chain by returning itself, while CommandExecutors consume one of the arguments
# to dispatch to a subcommand.
# eg of subcommands
#
# ./bin/my_app my_command help
# ./bin/my_app my_command subcommand
# class MyApp < Etna::Application
#   class MyCommand
#     include Etna::CommandExecutor
#
#     class Subcommand < Etna::Command
#       def execute
#       end
#     end
#   end
# end
#
module Etna
  # Provides the usage DSL and method that includes information about a CommandExecutor or a Command.
  # Commands or CommandExecutors can call 'usage' in their definition to provide a specific description
  # to be given about that command.  By default, the command name + desc method will be shown.
  module HasUsage
    module Dsl
      def usage(desc)
        define_method :usage do
          "  #{"%-45s" % command_name}#{desc}"
        end
      end
    end

    def command_name
      self.class.name.snake_case.split(/::/).last
    end

    def usage
      "  #{"%-45s" % command_name}#{desc}"
    end

    # By default, the description of a command maps the execute parameters into CLI descriptions,
    # where as CommandExecutors will display their subcommands.
    def desc
      if respond_to?(:execute)
        method(:execute).parameters.map do |type, name|
          name = "..." if name.nil?

          case type
          when :req
            "<#{name}>"
          when :opt
            "[<#{name}>]"
          when :rest
            "<#{name}>..."
          else
            raise "Invalid command execute argument specification, unsure how to format description."
          end
        end.join(' ')
      elsif respond_to?(:subcommands)
        '<command> <args>...'
      end
    end

    def self.included(cls)
      cls.extend(Dsl)
    end
  end

  # Include this module into class that dispatches to child CommandExecutors | Commands that must exist as
  # inner classes.  Note that non-root CommandExecutors must accept and set their @parent just like commands.
  module CommandExecutor
    # Program name is used to compose the usage display.  The top level executor will just assume the running
    # process's bin name, while other executors will use the 'command name' of the Command / Exectuor
    # (class name derived) and whatever program_name exists for its parent scope.
    def program_name
      if parent.nil? || !parent.respond_to?(:program_name)
        $PROGRAM_NAME
      else
        parent.program_name + " " + command_name
      end
    end

    attr_reader :parent
    def initialize(parent = nil)
      super()
      @parent = parent
    end

    def self.included(cls)
      cls.include(HasUsage)
      cls.const_set(:Help, Class.new(Etna::Command) do
        usage 'List this help'

        def execute
          self.parent.help
        end
      end) unless cls.const_defined?(:Help)
    end

    def help
      puts "usage: #{program_name} #{desc}"
      subcommands.each do |name, cmd|
        puts cmd.usage
      end
    end

    def find_command(cmd = 'help', *args)
      unless subcommands.include?(cmd)
        cmd = 'help'
        args = []
      end

      subcommands[cmd].find_command(*args)
    end

    def subcommands
      @subcommands ||= self.class.constants.sort.reduce({}) do |acc, n|
        acc.tap do
          c = self.class.const_get(n)
          next unless c.instance_methods.include?(:find_command)
          v = c.new(self)
          acc[v.command_name] = v
        end
      end
    end
  end

  class Command
    include HasUsage

    attr_reader :parent
    def initialize(parent = nil)
      @parent = parent
    end

    def self.parent_scope
      parts = self.name.split('::')
      parts.pop
      Kernel.const_get(parts.join('::'))
    end

    def find_command(*args)
      [self, args]
    end

    def fill_in_missing_params(args)
      req_params = method(:execute).parameters.select { |type, name| type == :req }
      args + (req_params[(args.length)..(req_params.length)] || []).map do |type, name|
        puts "#{name}?"
        STDIN.gets.chomp
      end
    end

    def completions
      method(:execute).parameters.map do |type, name|
        name = "..." if name.nil?
        if type == :req
          [completions_for(name)]
        else
          [[]]
        end
      end.inject([], &:+)
    end

    def completions_for(parameter)
      if parameter == 'env' || parameter == 'environment' || parameter =~ /_env/
        ['production', 'staging', 'development']
      else
        []
      end
    end


    # To be overridden during inheritance.
    def execute
      raise 'Command is not implemented'
    rescue => e
      Rollbar.error(e)
      raise
    end

    # To be overridden during inheritance, to e.g. connect to a database.
    # Should be called with super by inheriting method.
    def setup(config)
      Etna::Application.find(self.class).configure(config)
    rescue => e
      Rollbar.error(e)
      raise
    end
  end

  # application.rb instantiates this for the project scoping.
  # This generates a file, project-name.completion, which is sourced
  # in build.sh to provide autocompletion in that environment.
  class GenerateCompletionScript < Etna::Command
    def generate_for_scope(scope, lines, depth = 1)
      lines << "if [ \"${COMP_CWORD}\" == \"#{depth}\" ]; then"
      lines << "all_completion_names='#{scope.subcommands.keys.join(' ')}'"
      lines << 'COMPREPLY=($(compgen -W "$all_completion_names" "${COMP_WORDS[COMP_CWORD]}"))'
      lines << "else"
      scope.subcommands.each do |name, command|
        lines << "if [ \"${COMP_WORDS[#{depth}]}\" == '#{name}' ]; then"
        if command.class.included_modules.include?(CommandExecutor)
          generate_for_scope(command, lines, depth + 1)
        else
          completions = command.completions
          completions.each_with_index do |c, i|
            lines << "if [ \"${COMP_CWORD}\" == \"#{i + depth + 1}\" ]; then"
            if c.instance_of?(Array) && !c.empty?
              lines << "all_completion_names='#{c.join(' ')}'"
              lines << 'COMPREPLY=($(compgen -W "$all_completion_names" "${COMP_WORDS[COMP_CWORD]}"))'
            else
              lines << 'COMPREPLY=()'
              lines << "return"
            end
            lines << "fi"
          end
        end
        lines << "true"
        lines << "fi"
      end
      lines << "fi"
    end

    def execute
      name = parent.class.name.split('::').last.snake_case

      lines = []
      lines << "#!/usr/bin/env bash"
      lines << "_#{name}_completions() {"
      lines << "local all_completion_names=''"
      generate_for_scope(parent, lines)
      lines << "}"
      lines << "complete -o default -F _#{name}_completions #{name}"

      File.open("#{name}.completion", 'w') { |f| f.write(lines.join("\n") + "\n") }
    end
  end
end
