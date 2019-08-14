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

  class Schema < Etna::Command
    usage 'Show the current database schema.'

    def execute
      Metis.instance.db.tap do |db|
        db.extension(:schema_dumper)
        puts db.dump_schema_migration
      end
    end

    def setup(config)
      super
      Metis.instance.setup_db
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

  class Audit < Etna::Command
    usage '<project_name> Audit project directory'

    COLORS={
      red: 31,
      green: 32,
      yellow: 33
    }

    def colorize(text, color)
      "\e[#{COLORS[color]}m#{text}\e[0m"
    end

    def note(obj)
      @found ||= {}
      @found[obj.location] = obj
    end

    def audit_files(files)
      files.each do |file|
        note(file)
        puts "#{file.file_path} => #{colorize(file.location, ::File.exists?(file.location) ? :green : :red)}"
      end
    end

    def audit_folders(folders)
      folders.each do |folder|
        note(folder)
        puts "#{::File.join(folder.folder_path)} => #{colorize(folder.location, ::Dir.exists?(folder.location) ? :green : :red)}"
        audit_files(folder.files)
        audit_folders(folder.folders)
      end
    end

    def audit_dirs(dirs)
      dirs.each do |dir|
        next if @found[dir]
        folder_path = ::File.join(
          dir.sub(/^.*files/,'').split(/\//).map do |f|
            Metis::File.unsafe_file_name(f)
          end
        )
        puts "#{colorize(dir,:red)} => #{colorize(folder_path, :green)}"
        audit_blobs(
          Dir.glob("#{dir}/*").select{|l| !::File.directory?(l)}
        )
        audit_dirs(
          Dir.glob("#{dir}/*").select{|l| ::File.directory?(l)}
        )
      end
    end

    def audit_blobs(blobs)
      blobs.each do |blob|
        next if @found[blob]
        file_path = ::File.join(
          blob.sub(/^.*files/,'').split(/\//).map do |f|
            Metis::File.unsafe_file_name(f)
          end
        )
        puts "#{colorize(blob,:red)} => #{colorize(file_path, :green)}"
      end
    end

    def execute(project_name)
      buckets = Metis::Bucket.where(
        project_name: project_name,
      ).all

      buckets.each do |bucket|
        # database files
        audit_files(
          Metis::File.where(
            project_name: project_name, bucket: bucket, folder_id: nil
          ).all
        )

        audit_folders(
          Metis::Folder.where(
            project_name: project_name, bucket: bucket, folder_id: nil
          ).all
        )

        puts 'unaccounted:'
        audit_dirs(
          Dir.glob("#{bucket.location}/*").select{|l| ::File.directory?(l)}
        )

        audit_blobs(
          Dir.glob("#{bucket.location}/*").select{|l| !::File.directory?(l)}
        )
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end

  class Archive < Etna::Command
    usage 'Checksum and archive files.'

    def execute
      needs_hash = Metis::File.where(file_hash:nil).order(:updated_at).all[0..10]
      puts "Found #{needs_hash.count} files to be checksummed."
      needs_hash.each(&:compute_hash!)

      needs_archive = Metis::File.exclude(file_hash: nil).where(archive_id: nil).order(:updated_at).all[0..10]
      puts "Found #{needs_archive.count} files to be archived."
      needs_archive.each do |file|
        begin
          file.backup!
        rescue ArgumentError => e
          puts "Could not archive #{file.file_name}"
          next
        end
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end

  class Assimilate < Etna::Command
    usage '<project> <bucket> <path> <files> Add the indicated (non-tracked) files to Metis at the given path'

    def execute project_name, bucket_name, folder_path, *files
      folder_path = folder_path.sub(%r!^/!,'')

      bucket = Metis::Bucket.where(project_name: project_name, name: bucket_name).first

      if folder_path.empty?
        files.each do |file_path|
          Metis::Folder.assimilate(file_path, bucket)
        end
      else
        metis_folder = Metis::Folder.from_path(bucket, folder_path).last
        unless metis_folder
          puts "No such folder #{folder_path}"
          exit
        end
        files.each do |file|
          metis_folder.assimilate(file)
        end
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end
end
