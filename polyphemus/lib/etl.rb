require_relative 'helpers'

class Polyphemus
  class EtlExecutor
    def self.inherited(subclass)
      subclass.include(Etna::CommandExecutor)

      self.const_set(:Run, Class.new(Etna::Command) do
        usage 'runs the etl process until no more active processing is currently available.'

        def execute
          while true
            return unless self.parent_scope::EtlClass.new.run_once
            logger.info("Continuing to process...")
          end

          logger.info("No more work detected, resting")
        end
      end) unless self.const_defined?(:Run)

      self.const_set(:Reset, Class.new(Etna::Command) do
        usage 'resets the cursor for this etl, so that next processing starts from the beginning of time.'

        def execute
          logger.warn("Resetting etl cursors!")
          self.parent_scope::EtlClass.new.cursor_group.reset_all!
        end
      end) unless self.const_defined?(:Reset)
    end
  end

  class EtlCommand
    include Etna::CommandExecutor

    def self.command_name
      'etl'
    end

    def subcommands
      @subcommands ||= Etl.descendants.inject({}) do |acc, etl|
        acc.tap do
          acc[etl]
        end

        self.const_set(:"#{etl.name}", Class.new(EtlExecutor) do
          usage "manages the #{etl.name}"
          const_set(:EtlClass, etl)
        end)
      end
    end
  end

  # Abstract base class meant to process batches of data based on a cursor, scanner, and process.
  class Etl
    # Will use default client based on environment unless override is set in constructor
    include WithEtnaClients
    include WithLogger

    def self.inherited(subclass)
      descendants << subclass
    end

    def self.descendants
      @descendants ||= []
    end

    attr_reader :cursor_group

    def initialize(cursor_group:, scanner:)
      @scanner = scanner
      @cursor_group = cursor_group
    end

    # Returns true iff a batch was processed as part of this iteration.
    def run_once
      @cursor_group.with_next do |cursor|
        batch = @scanner.find_batch(cursor)

        if batch.empty? || process(cursor, batch) == :stop
          cursor.load_from_db
          return false
        end

        return cursor.save_to_db
      end

      false
    end

    # Subclasses should override with their processing loop.
    # IMPORTANT!  Keep in mind that because the process can fail part way through work, there is no guarantee of single
    # delivery and that process should be idempotent to the best of its ability.
    # process should return :stop when it wishes to abort processing the current batch, as may happen if the data
    # within the batch is not ready.  In this case, the cursor is not saved and a future run will again attempt to process
    # the batch.
    def process(cursor, batch) end

    def setup(config)
      super
      Polyphemus.instance.setup_db
    end
  end
end
