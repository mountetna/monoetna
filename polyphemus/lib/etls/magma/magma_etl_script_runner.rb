require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'

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

  def run(record, magma_client:, mode: :process)
    # Make the record available as an instance variable named based on the model_name.
    instance_variable_set(:"@#{model_name}", record)
    @record = record
    @magma_client = magma_client

    run_script

    if mode == :process
      process_outputs
    end
  end
end
