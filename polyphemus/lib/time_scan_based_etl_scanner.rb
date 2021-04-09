# Abstract base class for an etl that processes incoming data that can be sorted by updated at, but for which
# the precision of the updated_at is not strictly monotonic (multiple entries can belong to the same updated_at) and
# whose ordering can also change (but only in that entries can move forward in position in the ordering, never backwards).
require_relative 'etl_scanner'

class Polyphemus
  class TimeScanBasedEtlScanner < Polyphemus::EtlScanner
    def normalize_seen_ids(cursor)
      # Upgrade older versions of seen_ids.  This may duplicate work, but not a huge deal.
      unless cursor[:seen_ids].all? { |v| v.is_a?(Array) }
        cursor[:seen_ids] = []
      end
    end

    # Corresponds with the seen_ids and stability of a scan.
    # Should uniquely identify an update to a record.
    def result_scanned_position(r)
      [self.result_id(r), self.result_updated_at(r).to_i].freeze
    end

    def find_batch(cursor)
      expansion_number = 1
      cursor[:seen_ids] ||= []
      normalize_seen_ids(cursor)
      state = self.start_batch_state(cursor)
      last_scanned_positions = []
      last_updated_at = cursor.updated_at || Time.at(0)

      while true
        results = self.execute_batch_find(state, expansion_number)

        # The list of tuples identifying the records by both their id, and updated at,
        # to compare when we see a stable set of results.
        scanned_positions = results
              .map { |r| self.result_scanned_position(r) }
              .sort

        # Find the most updated item in the result set, and partition our results into those that
        # occur before that time, and those that do not
        last_updated_at = results.map { |r| self.result_updated_at(r) }.max unless results.empty?
        discrete_time_results, ongoing_time_set = results.partition do |result|
          !time_precision_overlaps?(self.result_updated_at(result), last_updated_at)
        end

        ongoing_time_set = ongoing_time_set.select { |v| !cursor[:seen_ids].include?(self.result_scanned_position(v)) }

        # After a previous scan, in which we did not yield a discrete dataset, if the last_scanned_positions
        # is stable with the next one (despite a larger expansion number), we assume that we are not going to
        # find a discrete time in the future, and thus commit to returning the ongoing_time_set, updating our seen_ids.
        if scanned_positions == last_scanned_positions
          cursor[:seen_ids] += ongoing_time_set.map { |r| self.result_scanned_position(r) }

          # Ensure that we have the most recent updated_at entry for any given seen id pair, as we can run into the
          # same id but under a new updated at.  We sort, and reverse, which should put duplicates in order of most recent,
          # then uniq by the id (first item)
          cursor[:seen_ids].sort!
          cursor[:seen_ids].reverse!
          cursor[:seen_ids].uniq!(&:first)

          return ongoing_time_set
        end

        last_scanned_positions = scanned_positions

        if discrete_time_results.empty?
          # Increase out paging size, see if we get new results and find a discrete time set.
          expansion_number += 1
        else
          cursor[:seen_ids] = discrete_time_results.map { |r| self.result_scanned_position(r) }
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