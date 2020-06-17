require 'date'
require 'logger'

class Archimedes
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Archimedes.instance.commands.each do |name,cmd|
        puts cmd.usage
      end
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Archimedes instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Archimedes.instance.setup_db
      Archimedes.instance.setup_magma
    end
  end
end
