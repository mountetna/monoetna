require 'shellwords'

class Vulcan

  class CommandBuilder
    def initialize
      @commands = []
    end

    def add(command, *args)
      @commands << [:command, command, *args]
      self
    end

    def add_raw(command)
      @commands << [:raw, command]
      self
    end

    def pipe_to(command, *args)
      @commands << [:pipe, command, *args]
      self
    end

    def redirect_to(file)
      @commands << [:redirect, file]
      self
    end

    def append_to(file)
      @commands << [:append, file]
      self
    end

    def background
      @commands << [:background]
      self
    end

    def to_s
      result = []
      
      @commands.each_with_index do |cmd, i|
        type = cmd[0]
        args = cmd[1..]
        
        case type
        when :command
          # Add && between commands unless it follows a pipe/redirect/append/background
          if i > 0 && ![:pipe, :redirect, :append, :background].include?(@commands[i-1][0])
            result << "&&"
          end
          result << Shellwords.join(args)
        when :raw
          if i > 0 && ![:pipe, :redirect, :append, :background].include?(@commands[i-1][0])
            result << "&&"
          end
          result << args[0]
        when :pipe
          result << "|"
          result << Shellwords.join(args)
        when :redirect
          result << ">"
          result << Shellwords.escape(args[0])
        when :append
          result << ">>"
          result << Shellwords.escape(args[0])
        when :background
          result << "&"
        end
      end
      
      result.join(' ')
    end
  end

end