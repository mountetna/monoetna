# Abstract base class for an etl
class Polyphemus
  class EtlScanner
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
  end
end
