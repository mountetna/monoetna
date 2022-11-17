
class Gnomon
  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
    end

    def setup(config)
    end
  end
  class Console < Etna::Command
    usage 'Open a console with a connected Gnomon instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Gnomon.instance.setup_db
      Gnomon.instance.setup_sequel
    end
  end
end

