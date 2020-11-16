require 'date'
require 'logger'

class Archimedes
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
