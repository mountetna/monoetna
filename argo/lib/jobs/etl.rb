
class ETLJob
    attr_reader :config
  
    def initialize(config = {})
      @config = config
    end
  
    # Pre-condition method to check if ETL should proceed
    #
    # @return [Boolean] True if conditions are met, false otherwise
    def pre
      raise NotImplementedError, "#{self.class} must implement the pre_condition method"
    end
  
    # Process method containing the main ETL logic
    def process
      raise NotImplementedError, "#{self.class} must implement the process method"
    end
  
    # Post-condition method to perform actions after processing
    def post
      raise NotImplementedError, "#{self.class} must implement the post method"
    end
  
    def execute
      run_etl
    end
  
    private
  
    # Orchestrates the ETL job by invoking pre, process, and post
    def run_etl
      if pre
        process
        post
        true
      else
        false
      end
    rescue StandardError => e
      raise e
    end
  end
  