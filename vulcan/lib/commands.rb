require 'date'
require 'logger'
require 'readline'

class Vulcan
  class Schema < Etna::Command
    usage 'Show the current database schema.'

    def execute
      Vulcan.instance.db.tap do |db|
        db.extension(:schema_dumper)
        puts db.dump_schema_migration
      end
    end

    def setup(config)
      super
      Vulcan.instance.setup_db
    end
  end

  class Workflow
    include Etna::CommandExecutor

    class Run < Etna::Command
      string_flags << '--session-key'
      string_flags << '--output-file'
      string_flags << '--input-file'

      def execute(project_name, workflow_name, input_file: nil, output_file: nil, session_key: SecureRandom.uuid.hex.to_s)
        @project_name = project_name
        @workflow_name = workflow_name
        @key = session_key
        @input_json = input_file&.yield_self do |f| JSON.parse(File.read(f)) end

        while true
          orchestration.run_until_done!(storage)

          find_all_output_ui_queries.each do |query_name, source, path|
            run_ui_query(query_name, :show_data, source, path)
          end

          found_needed_input = false
          while (input = find_next_primary_input)
            val = collect_primary_input(input.id, input.type)
            session.define_user_input([:primary_inputs, input.id], val)
            found_needed_input = true
          end

          while (ui_query = find_next_input_ui_query)
            query_name, source, input_files = ui_query

            val = run_ui_query(query_name, :query_json, source.join('/'), *input_files)
            session.define_user_input(source, val)
            found_needed_input = true
          end

          break unless found_needed_input
        end

        if output_file
          puts "Saving session to #{output_file}"
          File.open(output_file,"w") do |f|
            f.puts session.as_json.to_json
          end
        end
      end

      def run_ui_query(query_name, backup, *args)
        query_name = query_name.gsub(/-/, '_').gsub('.cwl', '').to_sym
        unless respond_to?(query_name)
          puts "Could not find #{query_name}, falling back"
          query_name = backup
        end

        send(query_name, *args)
      end


      def show_data(source, path)
        puts "** #{source.join('/')}:"
        puts ::File.read(path)
        puts
      end

      def query_json(name, *args)
        JSON.parse(Readline.readline("json value for #{name.inspect}: ", true))
      end

      def query_string(name, *args)
        Readline.readline("string value for #{name.inspect}: ", true)
      end

      def query_float(name, *args)
        while true
          val = Readline.readline("float value for #{name.inspect}: ", true)
          begin
            return Float(val)
          rescue
          end
        end
      end

      def query_int(name, *args)
        while true
          val = Readline.readline("int value for #{name.inspect}: ", true)
          begin
            return Integer(val)
          rescue
          end
        end
      end

      def collect_primary_input(id, type)
        case type
        when 'null'
          nil
        when 'string'
          query_string("primary input #{id}")
        when 'int'
          query_int("primary input #{id}")
        when 'long'
          query_float("primary input #{id}")
        when 'double'
          query_float("primary input #{id}")
        when 'float'
          query_float("primary input #{id}")
        else
          query_json("primary input #{id}")
        end
      end

      def find_next_primary_input
        workflow.inputs.each do |input|
          return input unless session.include?([:primary_inputs, input.id])
        end

        nil
      end

      def find_next_input_ui_query
        orchestration.next_runnable_build_targets(storage).each do |step_name, build_target|
          step = workflow.find_step(step_name)
          next if step.nil?

          query_name = step.ui_query_name
          next if query_name.nil?

          build_target.build_outputs.keys.each do |build_output_name|
            source = [step_name, build_output_name]
            return [query_name, source, build_target.input_files] unless session.include?(source)
          end
        end

        nil
      end

      def find_all_output_ui_queries
        [].tap do |result|
          orchestration.next_runnable_build_targets(storage).each do |step_name, build_target|
            step = workflow.find_step(step_name)

            if step
              query_name = step.ui_query_name
              next if query_name.nil?
              # "output" ui steps, other than primary outputs, are steps that lack any output
              next unless build_target.build_outputs.length == 0
            elsif step_name == :primary_outputs
              query_name = 'show_data'
            else
              next
            end

            if build_target.is_buildable?(storage) && build_target.is_built?(storage)
              build_target.input_files.each do |sf|
                result << [query_name, [step_name, sf.logical_name], sf.data_path(storage)]
              end
            end
          end
        end
      end

      def workflow
        session.workflow.tap do |workflow|
          if workflow.nil?
            raise "Workflow by name #{@workflow_name} does not exist!"
          end
        end
      end

      def session
        @session ||= Session.from_json(
          @input_json || {
            'project_name' => @project_name,
            'workflow_name' => @workflow_name,
            'key' => @key
          }
        )
      end

      def orchestration
        session.orchestration.tap do |orchestration|
          if orchestration.nil?
            raise Etna::NotFound.new("Workflow by the name #{session.workflow_name} could not be found.")
          end
        end
      end

      def storage
        @storage ||= Vulcan::Storage.new
      end

      def setup(config)
        super
        Vulcan.instance.setup_db
      end
    end
  end

  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Vulcan.instance.db

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
      Vulcan.instance.setup_db
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Vulcan instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Vulcan.instance.setup_db
    end
  end

  class CreateDb < Etna::Command
    usage '# create the initial database per config.yml'

    def execute
      if @no_db
        create_db if @no_db

        puts "Database is setup. Please run `bin/timur migrate #{@project_name}`."
      else
        puts "Database already exists."
      end
    end

    def create_db
      # Create the database only

      # We can't set environment vars for the executing shell,
      #   but we can set env vars for subshells that we open
      # https://stackoverflow.com/a/2660833
      fork do
        puts "Creating database #{@db_config[:database]}"
        ENV['PGPASSWORD'] = @db_config[:password]
        ENV['PGUSER'] = @db_config[:user]
        ENV['PGDB'] = @db_config[:database]
        %x{ createdb -w -U "$PGUSER" "$PGDB" }

        Vulcan.instance.setup_db
        exit 0
      end
      Process.wait
    end

    def setup(config)
      super
      @db_config = Vulcan.instance.config(:db)
      begin
        Vulcan.instance.setup_db
      rescue Sequel::DatabaseConnectionError
        @no_db = true
      end
    end
  end
end
