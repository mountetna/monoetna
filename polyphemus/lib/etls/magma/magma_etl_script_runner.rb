require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'
require_relative '../../data_processing/flow_jo_dsl'

class MagmaEtlScriptRunner < EtlScriptRunner
  include FlowJoDsl

  attr_reader :record, :magma_client
  def model_name
    if @script_name_parts.length < 3
      raise "Improperly formatted script name #{@file_path}.  Script name should start with '<project name>+<model name>+'."
    end

    @script_name_parts[1]
  end

  def process_outputs
    # TODO
  end

  def debug_outputs
    pp({ flow_jo_outputs: flow_jo_outputs })
  end

  def run(record, magma_client:, mode: :process)
    # Make the record available as an instance variable named based on the model_name.
    instance_variable_set(:"@#{model_name}", record)
    @record = record
    @magma_client = magma_client

    run_script

    if mode == :process
      process_outputs
    elsif mode == :debug
      debug_outputs
    end
  end
end
