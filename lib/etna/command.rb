require 'rollbar'

module Etna
  class Command
    class << self
      def usage(desc)
        define_method :usage do
          "  #{"%-45s" % name}#{desc}"
        end

        define_method :completions do
          []
        end
      end
    end

    def self.inherited(subclass)
      @descendants ||= []
      @descendants << subclass
    end

    def self.descendants
      @descendants ||= []
    end

    def usage
      "  #{"%-45s" % name}#{desc}"
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

    def desc
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
    end

    def name
      self.class.name.snake_case.split(/::/).last.to_sym
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

  class GenerateCompletionScript < Etna::Command
    def execute
      name = Etna::Application.instance.class.name.downcase
      commands = Etna::Application.instance.commands

      lines = []
      lines << "#!/usr/bin/env bash"
      lines << "_#{name}_completions() {"
      lines << "local all_completion_names=''"
      lines << 'if [ "${COMP_CWORD}" == "1" ]; then'
      lines << "all_completion_names='#{commands.map { |name, _| name }.join(' ')}'"
      lines << 'COMPREPLY=($(compgen -W "$all_completion_names" "${COMP_WORDS[COMP_CWORD]}"))'
      lines << "else"
      commands.each do |name, command|
        lines << "if [ \"${COMP_WORDS[1]}\" == '#{name}' ]; then"
        completions = command.completions
        completions.each_with_index do |c, i|
          lines << "if [ \"${COMP_CWORD}\" == \"#{i + 2}\" ]; then"
          if c.instance_of?(Array) && !c.empty?
            lines << "all_completion_names='#{c.join(' ')}'"
            lines << 'COMPREPLY=($(compgen -W "$all_completion_names" "${COMP_WORDS[COMP_CWORD]}"))'
          else
            lines << 'COMPREPLY=($(compgen -f "${COMP_WORDS[COMP_CWORD]}"))'
            lines << "return"
          end
          lines << "fi"
        end
        lines << "true"
        lines << "fi"
      end
      lines << 'fi'
      lines << "}"
      lines << "complete -F _#{name}_completions #{name}"

      File.open("#{name}.completion", 'w') { |f| f.write(lines.join("\n") + "\n") }
    end
  end

  def completions_for(parameter)
    case parameter
    when "output_file"
    end
  end
end
