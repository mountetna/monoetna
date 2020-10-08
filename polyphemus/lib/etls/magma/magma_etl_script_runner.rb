require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'

class MagmaEtlScriptRunner < EtlScriptRunner
  include XMLDsl
  include MagmaDsl
  include FlowJoDsl

  def model_name
    if @script_name_parts.length < 3
      raise "Improperly formatted script name #{@file_path}.  Script name should start with '<project name>_<model name>_'."
    end

    @script_name_parts[1]
  end

  def process_outputs
    # TODO
  end

  def run(record)
    # Make the record available as an instance variable named based on the model_name.
    instance_variable_set(:"@#{model_name}", record)
    run_script

    process_outputs
  end
end
