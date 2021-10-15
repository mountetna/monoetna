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

      subclass.const_set(:Run, Class.new(Etna::Command) do
        usage 'runs the etl process until no more active processing is currently available.'
        include WithLogger

        def etl
          @etl ||= self.parent.class::EtlClass.new
        end

        def execute
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

      subclass.const_set(:Reset, Class.new(Etna::Command) do
        usage 'resets the cursor for this etl, so that next processing starts from the beginning of time.'
        include WithLogger

        def etl
          @etl ||= self.parent.class::EtlClass.new
        end

        def execute
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

    def self.inherited(subclass)
      @descendants ||= []
      @descendants << subclass
    end

    def self.descendants
      @descendants ||= []
      @descendants + @descendants.map(&:descendants).inject([], &:+)
    end

    attr_reader :cursor_group

    def initialize(cursor_group:, scanner:)
      @scanner = scanner
      @cursor_group = cursor_group
    end

    # Returns true iff a batch was processed as part of this iteration.
    def run_once
      logger.info("Starting loop")
      @cursor_group.with_next do |cursor|
        logger.info("Selecting cursor for #{cursor.name}")
        cursor.load_from_db
        logger.info("Finding batch...")
        batch = @scanner.find_batch(cursor)

        logger.info("Attempting to process")
        if batch.empty? || process(cursor, batch) == :stop
          logger.info("No more work found, stopping.")
          next false
        end

        logger.info("Saving cursor to database")
        cursor.save_to_db
        next true
      end
    end

    # Subclasses should override with their processing loop.
    # IMPORTANT!  Keep in mind that because the process can fail part way through work, there is no guarantee of single
    # delivery and that process should be idempotent to the best of its ability.
    # process should return :stop when it wishes to abort processing the current batch, as may happen if the data
    # within the batch is not ready.  In this case, the cursor is not saved and a future run will again attempt to process
    # the batch.
    def process(cursor, batch) end
  end
end
