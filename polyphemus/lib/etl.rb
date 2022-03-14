require_relative "helpers"

class Polyphemus
  class EtlExecutor
    def self.ensure_for_etl(etl)
      return nil unless etl.descendants.empty? # Ignore base classes.
      name = etl.name.split("::").last
      EtlCommand.const_set(:"#{name}", Class.new(EtlExecutor) do
        usage "manages the #{name}"
        const_set(:EtlClass, etl)
      end) unless EtlCommand.constants.include?(:"#{name}")
    end

    def self.inherited(subclass)
      subclass.include(Etna::CommandExecutor)
      subclass.boolean_flags << '--from-environment'

      subclass.const_set(:Run, Class.new(Etna::Command) do
        usage 'runs the etl process until no more active processing is currently available.'
        include WithLogger

        def etl
          if @from_environment
            @etl ||= self.parent.class::EtlClass.from_hash(Etna::EnvironmentVariables.load_from_env('ETL') { |path, value| [path, YAML.load(value)] }, true)
          else
            @etl ||= self.parent.class::EtlClass.new
          end
        end

        def execute(from_environment: false)
          @from_environment = from_environment

          while true
            break unless etl.run_once
            logger.info("Continuing to process...")
          end

          logger.info("No more work detected, resting")
        end

        def setup(config)
          super
          Polyphemus.instance.setup_db
          Polyphemus.instance.setup_logger
          Polyphemus.instance.setup_sequel
          Polyphemus.instance.setup_ssh
        end
      end) unless self.const_defined?(:Run)

      subclass.const_set(:FindBatch, Class.new(Etna::Command) do
        usage 'selects a single batch of work for this etl.'
        include WithLogger

        # Used in testing.  Not the best, but the commands are singletons which are 'instantiated' from the execute
        # call, but in our tests we access the .etl before execute is called.
        def enable_from_environment
          @from_environment = true
        end

        def etl
          if @from_environment
            @etl ||= self.parent.class::EtlClass.from_hash(
              Etna::EnvironmentVariables.load_from_env('ETL') { |path, value| [path, YAML.load(value)] }, true)
          else
            @etl ||= self.parent.class::EtlClass.new
          end
        end

        def execute(from_environment: false)
          @from_environment = from_environment

          batches = []

          while true
            break unless etl.run_once do |cursor, batch|
              logger.info("Consuming #{batch.length} items from #{cursor.updated_at}...")
              batches.push(*batch)
            end
            logger.info("Continuing to process...")
          end

          logger.info("No more batches found, outputting result")
          dump_result(etl.serialize_batch(batches))
        end

        def dump_result(result)
          puts(result.to_json)
        end

        def setup(config)
          super
          Polyphemus.instance.setup_db
          Polyphemus.instance.setup_logger
          Polyphemus.instance.setup_sequel
          Polyphemus.instance.setup_ssh
        end
      end) unless self.const_defined?(:FindBatch)

      subclass.const_set(:Reset, Class.new(Etna::Command) do
        usage 'resets the cursor for this etl, so that next processing starts from the beginning of time.'
        include WithLogger

        def etl
          @etl ||= self.parent.class::EtlClass.new
        end

        def execute(from_environment: false)
          if from_environment
            raise "Cannot reset etl cursors with --from-environment set!"
          end

          logger.warn("Resetting etl cursors!")
          etl.cursor_group.reset_all!
        end

        def setup(config)
          super
          Polyphemus.instance.setup_db
          Polyphemus.instance.setup_logger
        end
      end) unless self.const_defined?(:Reset)
    end
  end

  class EtlCommand
    include Etna::CommandExecutor

    def command_name
      "etl"
    end

    def subcommands
      Etl.descendants.each do |etl|
        EtlExecutor.ensure_for_etl(etl)
      end

      super
    end
  end

  class EtlError < Exception
  end

  # Abstract base class meant to process batches of data based on a cursor, scanner, and process.
  class Etl
    # Will use default client based on environment unless override is set in constructor
    include WithEtnaClients
    include WithLogger
    extend ::Etna::Injection::FromHash

    attr_reader :scanner, :cursor_group

    def self.inherited(subclass)
      @descendants ||= []
      @descendants << subclass
    end

    def self.descendants
      @descendants ||= []
      @descendants + @descendants.map(&:descendants).inject([], &:+)
    end

    def initialize(cursors:, scanner:)
      @cursor_group = EtlCursorGroup.new(cursors)
      @scanner = scanner
    end

    def cursors_from_pairs(pairs:, pair_keys:, cls:, cursor_env:)
      pairs.map do |pair|
        cls.from_hash({"job_name" => self.class.name}.update(pair_keys.zip(pair).to_h).update(cursor_env), true).tap do |cursor|
          cursor.load_from_db unless cursor.from_env?
        end
      end
    end

    # To support processing batches adding as inputs
    def serialize_batch(batch)
      batch.to_json
    end

    def deserialize_batch(string)
      JSON.parse(string, symbolize_names: true)
    end

    # Subclasses should override with their processing loop.
    # IMPORTANT!  Keep in mind that because the process can fail part way through work, there is no guarantee of single
    # delivery and that process should be idempotent to the best of its ability.
    # process should return :stop when it wishes to abort processing the current batch, as may happen if the data
    # within the batch is not ready.  In this case, the cursor is not saved and a future run will again attempt to process
    # the batch.
    def process(cursor, batch) end

    def find_batch(cursor)
      scanner.find_batch(cursor)
    end

    # Returns true iff a batch was processed as part of this iteration.
    def run_once(&processor)
      processor = Proc.new { |cursor, batch| process(cursor, batch) } unless block_given?

      logger.info("Starting loop")
      cursor_group.with_next do |cursor|
        logger.info("Selecting cursor for #{cursor.name}")
        cursor.load_from_db unless cursor.from_env?
        logger.info("Finding batch...")
        batch = find_batch(cursor)

        logger.info("Attempting to process")
        if batch.empty? || processor.call(cursor, batch) == :stop
          logger.info("No more work found, stopping.")
          next false
        end

        cursor.save_to_db unless cursor.from_env?
        next true
      end
    end
  end
end
