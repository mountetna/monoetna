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
  module CommandOrExecutor
    module Dsl
      def usage(desc)
        define_method :usage do
          "  #{"%-45s" % command_name}#{desc}"
        end
      end

      def boolean_flags
        @boolean_flags ||= []
      end

      def string_flags
        @string_flags ||= []
      end

      def multi_flags
        @multi_flags ||= []
      end
    end

    def flag_as_parameter(flag)
      flag.gsub('--', '').gsub('-', '_')
    end

    def parse_flags(*args)
      new_args = []
      flags = {}
      found_non_flag = false

      until args.empty?
        next_arg = args.shift

        unless next_arg.start_with? '--'
          new_args << next_arg
          found_non_flag = true
          next
        end

        arg_name = flag_as_parameter(next_arg).to_sym

        if self.class.boolean_flags.include?(next_arg)
          flags[arg_name] = true
        elsif self.class.string_flags.include?(next_arg)
          if args.empty?
            raise "flag #{next_arg} requires an argument"
          else
            flags[arg_name] = args.shift
          end
        elsif self.class.multi_flags.include?(next_arg)
          if args.empty?
            raise "flag #{next_arg} requires an argument"
          else
            (flags[arg_name] ||= []) << args.shift
          end
        elsif !found_non_flag
          raise "#{program_name} does not recognize flag #{next_arg}"
        else
          new_args << next_arg
        end
      end

      [flags, new_args]
    end

    def completions_for(parameter)
      if parameter == 'env' || parameter == 'environment' || parameter =~ /_env/
        ['production', 'staging', 'development']
      else
        ["__#{parameter}__"]
      end
    end

    def command_name
      self.class.name.snake_case.split(/::/).last
    end

    def usage
      "  #{"%-45s" % command_name}#{desc}"
    end

    def flag_argspec
      @argspec ||= []
    end

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
          when :keyrest
            "[flags...]"
          when :key
            flag = "--#{name.to_s.gsub('_', '-')}"
            if self.class.boolean_flags.include?(flag)
              "[#{flag}]"
            else
              "[#{flag} <#{name}>]"
            end
          else
            raise "Invalid command execute argument specification #{type}, unsure how to format description."
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
    attr_reader :parent

    def initialize(parent = nil)
      super()
      @parent = parent
    end

    def self.included(cls)
      cls.include(CommandOrExecutor)
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

    def find_command(*args, **kwds)
      flags, args = parse_flags(*args)
      dispatch_to_subcommand(*args, **(kwds.update(flags)))
    end

    def dispatch_to_subcommand(cmd = 'help', *args, **kwds)
      unless subcommands.include?(cmd)
        cmd = 'help'
        args = []
      end

      subcommands[cmd].find_command(*args, **kwds)
    end

    def subcommands
      @subcommands ||= self.class.constants.sort.reduce({}) do |acc, n|
        acc.tap do
          c = self.class.const_get(n)
          next unless c.respond_to?(:instance_methods)
          next unless c.instance_methods.include?(:find_command)
          v = c.new(self)
          acc[v.command_name] = v
        end
      end
    end

    def all_subcommands
      subcommands.values + (subcommands.values.map { |s| s.respond_to?(:all_subcommands) ? s.all_subcommands : [] }.flatten)
    end
  end

  class Command
    include CommandOrExecutor

    attr_reader :parent

    def initialize(parent = nil)
      @parent = parent
    end

    def self.parent_scope
      parts = self.name.split('::')
      parts.pop
      Kernel.const_get(parts.join('::'))
    end

    def find_command(*args, **kwds)
      flags, args = parse_flags(*args)
      [self, args, kwds.update(flags)]
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
        if type == :req || type == :opt
          [completions_for(name)]
        else
          []
        end
      end.inject([], &:+)
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
end
