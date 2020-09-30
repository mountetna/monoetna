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
      expansion_number = 1
      cursor[:seen_ids] ||= []
      state = self.start_batch_state(cursor)
      last_batch_ids = []
      last_updated_at = Time.at(0)

      while true
        results = self.execute_batch_find(state, expansion_number)
        prefilter_ids = results
              .map { |r| self.result_id(r) }
              .sort

        results = results
            .select { |r| !cursor[:seen_ids].include?(self.result_id(r)) }

        ids = results
            .map { |r| self.result_id(r) }
            .sort

        # No progress made
        if prefilter_ids == last_batch_ids
          cursor[:seen_ids] += ids
          cursor[:seen_ids].uniq!
          puts "Oh no! #{last_batch_ids.length} #{results.length}"
          return results
        end
        last_batch_ids = prefilter_ids

        last_updated_at = results.map { |r| self.result_updated_at(r) }.max unless results.empty?
        discrete_time_results = results.select do |result|
          !time_precision_overlaps?(self.result_updated_at(result), last_updated_at)
        end

        if discrete_time_results.empty?
          expansion_number += 1
        else
          cursor[:seen_ids] = discrete_time_results.map { |r| self.result_id(r) }
          cursor.updated_at = discrete_time_results.map { |r| self.result_updated_at(r) }.max
          puts "awww yeaaa #{cursor.updated_at}"
          return discrete_time_results
        end
      end
    end

    def time_precision_overlaps?(a_time, b_time)
      a_time.to_i == b_time.to_i
    end
  end
end