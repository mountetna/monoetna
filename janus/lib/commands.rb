class Janus
  class Console < Etna::Command
    usage 'Open a console with a connected magma instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      Janus.instance.configure(config)
      Janus.instance.setup_db
    end
  end

  class AddUser < Etna::Command
    def execute(email, name, password=nil)
      user = User.find_or_create(email: email) do |user|
        user.name = name.strip
      end

      user.tap do |user|
        user.name = name.strip
        if password
          user.pass_hash = Janus.instance.sign.hash_password(password)
        end
        user.save
      end
    end

    def setup(config)
      Janus.instance.configure(config)
      Janus.instance.setup_db
    end
  end

  class AddUserKey < Etna::Command
    usage '<email> <pem_file>'
    def execute(email, pem_file)
      user = User[email: email]
      if !user
        puts 'User not found.'
        exit
      end

      if !File.exists?(pem_file)
        puts 'No such key file'
        exit
      end

      pem = File.read(pem_file)
      begin
        key = OpenSSL::PKey::RSA.new(pem)
      rescue
        puts 'Could not parse key file!'
        exit
      end

      user.public_key = pem
      user.save
    end

    def setup(config)
      Janus.instance.configure(config)
      Janus.instance.setup_db
    end
  end

  class AddProject < Etna::Command
    usage '<project_name> <project_name_full>'
    def execute(project_name, project_name_full)
      attributes = { project_name: project_name }
      project = Project.find(attributes) || Project.new(attributes)
      project.project_name_full = project_name_full
      project.save
    end

    def setup(config)
      Janus.instance.configure(config)
      Janus.instance.setup_db
    end
  end

  class Permit < Etna::Command
    boolean_flags << '--privileged'

    def execute(email, project_name, role, privileged: false)
      user = User[email: email]
      project = Project[project_name: project_name]
      if !user
        puts 'User not found.'
        exit
      end

      if !project
        puts 'Project not found.'
        exit
      end

      attributes = { project: project, user: user, privileged: privileged }

      perm = Permission.find(attributes) || Permission.new(attributes)
      perm.role = role
      perm.save
    end

    def setup(config)
      Janus.instance.configure(config)
      Janus.instance.setup_db
    end
  end

  class Migrate < Etna::Command
    usage '[<version number>] # blank to migrate to the latest'
    string_flags << '--version'
    
    def execute(version: nil)
      Sequel.extension(:migration)
      db = Janus.instance.db

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
      Janus.instance.setup_db(false)
    end
  end

  class GenerateKeyPair < Etna::Command
    usage '<key_size> # Generate a private/public key pair in PEM format'

    def execute(key_size)
      key_size = key_size.to_i

      if key_size < 1024
        puts 'Your key size is too small'
        exit
      end

      if Math.log(key_size,2) != Math.log(key_size,2).to_i
        puts 'Key size must be a power of 2'
        exit
      end

      private_key = Janus.instance.sign.generate_private_key(key_size)

      puts 'Private key:'
      puts private_key

      puts

      puts 'Public key:'
      puts private_key.public_key
    end
  end
end
