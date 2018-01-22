class Metis
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Metis.instance.commands.each do |name,cmd|
        puts cmd.usage
      end
    end

    def setup(config)
      Metis.instance.configure(config)
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected magma instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      Metis.instance.configure(config)
      Metis.instance.load_models
    end
  end
  
  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    
    def execute(version=nil)
      Sequel.extension(:migration)
      db = Metis.instance.db

      if version
        puts "Migrating to version #{version}"
        Sequel::Migrator.run(db, 'db/migrations', target: version.to_i)
      else
        puts 'Migrating to latest'
        Sequel::Migrator.run(db, 'db/migrations')
      end
    end

    def setup(config)
      super
      Metis.instance.setup_db
    end
  end
end
