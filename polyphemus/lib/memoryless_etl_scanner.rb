# Abstract base class for an etl that processes incoming data and accepts it all.
# Useful in scenarios where no built-in scanning interface exists or is required, i.e.
#   rsync or SFTP.
require_relative "etl_scanner"

class Polyphemus
  class MemorylessEtlScanner < Polyphemus::EtlScanner
    def find_batch(cursor)
      cursor[:seen_ids] ||= []
      self.execute_batch_find
    end
  end
end
