# Abstract base class for an etl that processes incoming data that can be sorted by updated at, but for which
# the precision of the updated_at is not strictly monotonic (multiple entries can belong to the same updated_at) and
# whose ordering can also change (but only in that entries can move forward in position in the ordering, never backwards).
class Polyphemus
  class TimeScanBasedEtlScanner
    def initialize
    end

    def start_batch_state(*args, &block)
      if block.nil?
        @start_batch_state.call(*args)
      else
        @start_batch_state = block
        self
      end
    end

    def execute_batch_find(*args, &block)
      if block.nil?
        @execute_batch_find.call(*args)
      else
        @execute_batch_find = block
        self
      end
    end

    def result_updated_at(*args, &block)
      if block.nil?
        @result_updated_at.call(*args)
      else
        @result_updated_at = block
        self
      end
    end

    def result_id(*args, &block)
      if block.nil?
        @result_id.call(*args)
      else
        @result_id = block
        self
      end
    end

    def find_batch(cursor)
      offset = 0
      cursor[:seen_ids] ||= []
      state = self.start_batch_state(cursor)
      non_discrete_batch = []
      last_updated_at = Time.at(0)

      while true
        results = self.execute_batch_find(state, offset)
        last_updated_at_results = results.map { |r| self.result_updated_at(r) }.max
        last_updated_at = [last_updated_at_results, last_updated_at].max unless last_updated_at_results.nil?
        non_discrete_batch.push(*results.select { |r| !cursor[:seen_ids].include?(self.result_id(r)) })

        if results.empty?
          return non_discrete_batch
        end

        discrete_time_results = non_discrete_batch.select do |result|
          !time_precision_overlaps?(self.result_updated_at(result), last_updated_at)
        end

        if discrete_time_results.empty?
          offset += results.length
          cursor[:seen_ids] += results.map { |r| self.result_id(r) }
        else
          cursor[:seen_ids] = discrete_time_results.map { |r| self.result_id(r) }
          cursor.updated_at = discrete_time_results.map { |r| self.result_updated_at(r) }.max
          return discrete_time_results
        end
      end
    end

    def time_precision_overlaps?(a_time, b_time)
      a_time.to_i == b_time.to_i
    end
  end
end