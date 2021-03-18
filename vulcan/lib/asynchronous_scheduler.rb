class Vulcan
  class AsynchronousScheduler
    def initialize(orchestration:)
      @orchestration = orchestration
      @running = {}
      @running_semaphore = Mutex.new
      @orchestration_semaphore = Mutex.new
    end

    def workflow
      @orchestration.workflow
    end

    def session
      @orchestration.session
    end

    def labels
      @labels ||= {
          hash: 'vulcan-hash',
          token: 'vulcan-token',
      }
    end

    def join_all
      threads = @running_semaphore.synchronize do
        @running.values
      end

      threads.each do |t|
        t.join
      end
    end

    def running?(build_target)
      @running_semaphore.synchronize do
        if @running.include?(build_target.cell_hash)
          return {status: 'running'}
        end
      end
    end

    #
    # def status(storage:, build_target:)
    #   hash = build_target.cell_hash
    #
    #   if build_target.is_built?(storage)
    #     return {status: 'complete'}
    #   end
    #
    #   @running_semaphore.synchronize do
    #     if @running.include?(hash)
    #       return {status: 'running'}
    #     end
    #   end
    #
    #   if (err = WorkflowError.find_error(hash: hash))
    #     return {status: 'error', error: err.message}
    #   end
    #
    #   return {status: 'pending'}
    # end

    # Schedules a thread to build the given build_target, thn runs it via orchestration, cleaning up
    # as it goes.
    def schedule!(build_target:, token:, storage:)
      hash = build_target.cell_hash
      # Ensure that any previous errors are cleared as we start a new process
      s = @running_semaphore

      s.synchronize do
        puts "inside schedule, checking included... #{hash}"

        if @running.include?(hash)
          return 0
        end

        puts "not included #{hash}"

        WorkflowError.find_error(hash: hash)&.clear!

        @running[hash] = Thread.new do
          begin
            puts "trying to run for hash #{hash}"
            @orchestration.run!(storage: storage, build_target: build_target, token: token)
            puts "Finished running, now scheduling more..."
            schedule_more!(token: token, storage: storage)
          rescue => e
            puts "Error #{e.to_s}"
            WorkflowError.mark_error!(hash: hash, message: e.to_s)
          ensure
            s.synchronize do
              puts "Removing #{hash} from running"
              @running.delete(hash)
            end
          end
        end
      end

      1
    end

    # Finds any runnable jobs, clears errors associated with their hashes, and starts processes for them.
    def schedule_more!(token:, storage:)
      @orchestration_semaphore.synchronize do
        puts "scheduling more..."
        runnables = @orchestration.next_runnable_build_targets(storage)

        started = 0

        runnables.each do |build_target|
          started += schedule!(build_target: build_target, token: token, storage: storage)
        end

        started
      end
    end
  end
end
