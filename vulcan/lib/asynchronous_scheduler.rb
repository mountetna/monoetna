class Vulcan
  class AsynchronousScheduler
    def initialize
      @running = {}
      @running_semaphore = Mutex.new
      @orchestration_semaphore = Mutex.new

      # Used in testing to control when jobs actually start to run, in order to allow better consistent test behavior.
      # Not exposed or used anywhere else.
      @_test_gating_semaphore = Mutex.new
    end

    def labels
      @labels ||= {
          hash: 'vulcan-hash',
          token: 'vulcan-token',
      }
    end

    def join_scheduled
      threads = @running_semaphore.synchronize do
        @running.values
      end

      threads.each do |t|
        t.join
      end

      threads.length
    end

    def join_all
      while join_scheduled > 0
      end
    end

    def running?(build_target)
      @running_semaphore.synchronize do
        if @running.include?(build_target.cell_hash)
          return {status: 'running'}
        end
      end
    end

    def status(storage:, build_target:, step: nil)
      hash = build_target.cell_hash

      requires_inputs = step&.ui_output_name && !build_target.is_buildable?(storage)

      if !requires_inputs && build_target.is_built?(storage)
        return {status: 'complete'}
      end

      @running_semaphore.synchronize do
        if @running.include?(hash)
          return {status: 'running'}
        end
      end

      if (err = WorkflowError.find_error(hash: hash))
        return {status: 'error', error: err.message}
      end

      return {status: 'pending'}
    end

    # Schedules a thread to build the given build_target, thn runs it via orchestration, cleaning up
    # as it goes.
    def schedule!(orchestration:, build_target:, token:, storage:)
      hash = build_target.cell_hash
      # Ensure that any previous errors are cleared as we start a new process
      s = @running_semaphore

      s.synchronize do
        Vulcan.instance.logger.info "inside schedule, checking included... #{hash}"

        if @running.include?(hash)
          return false
        end

        Vulcan.instance.logger.info "not included #{hash}"

        WorkflowError.find_error(hash: hash)&.clear!

        @running[hash] = Thread.new do
          begin
            @_test_gating_semaphore.synchronize do
              Vulcan.instance.logger.info "trying to run for hash #{hash}"
              orchestration.run!(storage: storage, build_target: build_target, token: token)
              Vulcan.instance.logger.info "Finished running, now scheduling more..."
              schedule_more!(orchestration: orchestration, token: token, storage: storage)
            end
          rescue => e
            Vulcan.instance.logger.error "Error #{e.to_s}"
            WorkflowError.mark_error!(hash: hash, message: e.to_s)
          ensure
            s.synchronize do
              Vulcan.instance.logger.info "Removing #{hash} from running"
              @running.delete(hash)
            end
          end
        end
      end

      true
    end

    # Finds any runnable jobs, clears errors associated with their hashes, and starts processes for them.
    def schedule_more!(token:, storage:, orchestration:)
      @orchestration_semaphore.synchronize do
        Vulcan.instance.logger.info("scheduling more...")
        runnables = orchestration.next_runnable_build_targets(storage)

        started = []

        runnables.each do |build_target|
          started << build_target.cell_hash if schedule!(orchestration: orchestration, build_target: build_target, token: token, storage: storage)
        end

        started
      end
    end
  end
end
