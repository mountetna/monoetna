require_relative "etl_scanner"

# Abstract base class that uses hash scan to determine what records to update
class Polyphemus
  class HashScanBasedEtlScanner < Polyphemus::EtlScanner
    def result_file_hashes(*args, &block)
      if block.nil?
        @result_file_hash.call(*args)
      else
        @result_file_hash = block
        self
      end
    end

    # Corresponds with the seen_ids and stability of a scan.
    # Should uniquely identify an update to a record.
    def result_scanned_position(r)
      [self.result_id(r),
       self.result_file_hashes(r)].freeze
    end

    def find_batch(cursor)
      expansion_number = 1
      cursor[:seen_ids] ||= []
      state = self.start_batch_state(cursor)
      last_scanned_positions = []
      last_updated_at = cursor.updated_at || Time.at(0)

      while true
        results = self.execute_batch_find(state, expansion_number)

        # The list of tuples identifying the records by both their id, and file hashes,
        # to compare when we see a stable set of results.
        scanned_positions = results
          .map { |r| self.result_scanned_position(r) }
          .sort

        # Partition our results into those that
        # have a changed hash and those who do not.
        changed_set, unchanged_set = results.partition do |result|
          file_hashes_changed?(cursor, result)
        end

        # After a previous scan, in which we did not yield a changed dataset, if the last_scanned_positions
        # is stable with the next one (despite a larger expansion number), we assume that we are not going to
        # find a changed set in the future, and thus commit to returning the changed_set, updating our seen_ids.
        if scanned_positions == last_scanned_positions
          update_seen_ids(cursor, changed_set)

          return changed_set
        end

        last_scanned_positions = scanned_positions

        if changed_set.empty?
          # Increase out paging size, see if we get new results and find a changed set.
          expansion_number += 1
        else
          # Persist all scanned records here with their last known
          #  file hashes, so we can compare them against future changes.
          update_seen_ids(cursor, changed_set)
          cursor.updated_at = changed_set.map { |r| self.result_updated_at(r) }.max
          return changed_set
        end
      end
    end

    def update_seen_ids(cursor, changed_set)
      # Drop any previously seen ids that are in the new changed set,
      #   so we will have the most recent hash entry for any given seen id pair.
      changed_set_ids = changed_set.map { |r| self.result_id(r) }
      cursor[:seen_ids] = cursor[:seen_ids].select do |seen_id|
        !changed_set_ids.include?(seen_id.first)
      end
      cursor[:seen_ids] += changed_set.map { |r| self.result_scanned_position(r) }
    end

    def file_hashes_changed?(cursor, result)
      # File hash changed if the result scanned is not present
      !cursor[:seen_ids].include?(self.result_scanned_position(result))
    end
  end
end
