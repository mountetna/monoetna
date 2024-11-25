
class ETLJob
    attr_reader :config
    
  
    def initialize(config = {}, secrets = {})
      @config = config
      @secrets = secrets
    end
    
    def run_id
      @run_id ||= ENV['ARGO_WORKFLOW_ID']
    end

    def workflow_table
      @workflow_table ||= ENV['WORKFLOW_TABLE']
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