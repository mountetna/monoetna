require_relative 'helpers'
class Polyphemus
  # Abstract base class meant to process batches of data based on a cursor, scanner, and process.
  class Etl < Etna::Command
    # Will use default client based on environment unless override is set in constructor
    include WithEtnaClients
    include WithLogger

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
    def process(cursor, batch)
    end

    def execute(subcommand)
      if subcommand == 'reset'
        logger.warn("Resetting etl cursors!")
        @cursor_group.reset_all!
      elsif subcommand == 'run'
        while true
          return unless run_once
          logger.info("Continuing to process...")
        end

        logger.info("No more work detected, resting")
      else
        raise "Unknown subcommand #{subcommand}, did you mean #{completions_for('subcommand').join(' | ')}?"
      end
    end

    def completions_for(parameter)
      if parameter == 'subcommand'
        ['reset', 'run']
      else
        []
      end
    end

    def setup(config)
      super
      Polyphemus.instance.setup_db
    end
  end
end
