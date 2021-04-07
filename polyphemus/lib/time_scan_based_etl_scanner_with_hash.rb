require_relative "time_scan_based_etl_scanner"

# Abstract base class that uses time scan to determine what records to update
#   as well as stores a file reference's file hash
class Polyphemus
  class TimeScanBasedEtlScannerWithHash < Polyphemus::TimeScanBasedEtlScanner
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
       self.result_updated_at(r).to_i,
       self.result_file_hashes(r)].freeze
    end
  end
end
