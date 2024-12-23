class Polyphemus
  class ETLJob
    attr_reader :config
  
    def initialize(config = {}, runtime_config = {}) 
      @config = config # config is the entire record for the table configs
      @runtime_config = runtime_config # runtime_config is the runtime_config.config json column
    end
    
    def run_id
      @run_id ||= ENV['ARGO_WORKFLOW_ID']
    end

    # Pre-condition method to check if ETL should proceed
    #
    # @return [Boolean] True if conditions are met, false otherwise
    def pre(context)
      raise NotImplementedError, "#{self.class} must implement the pre_condition method"
    end
  
    # Process method containing the main ETL logic
    def process(context)
      raise NotImplementedError, "#{self.class} must implement the process method"
    end
  
    # Post-condition method to perform actions after processing
    def post(context)
      raise NotImplementedError, "#{self.class} must implement the post method"
    end
  
    def execute
      run_etl
    end
  
    private
  
    # Orchestrates the ETL job by invoking pre, process, and post
    def run_etl
      context = {}
      if pre(context)
        process(context)
        post(context)
      end
    rescue StandardError => e
      raise e
    end
  end
end
