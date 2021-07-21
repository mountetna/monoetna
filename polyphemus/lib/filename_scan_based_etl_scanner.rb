# Abstract base class for an etl that processes incoming data and scans
#   based on file names.
# Useful in scenarios where no time or hash-based scanning is available, i.e.
#   rsync.
require_relative "etl_scanner"

class Polyphemus
  class FilenameScanBasedEtlScanner < Polyphemus::EtlScanner
    # Corresponds with the seen_ids and stability of a scan.
    # Should uniquely identify an update to a record.
    def result_scanned_position(r)
      [self.result_id(r), self.result_updated_at(r).to_s].freeze
    end

    def find_batch(cursor)
      cursor[:seen_ids] ||= []
      results = self.execute_batch_find

      # The list of tuples identifying the records by both their id, and updated at,
      # to compare when we see a stable set of results.
      scanned_positions = results
        .map { |r| self.result_scanned_position(r) }
        .sort
        .reverse
        .uniq(&:first)

      return [] if scanned_positions == cursor[:seen_ids]

      cursor[:seen_ids] = scanned_positions

      results
    end
  end
end
