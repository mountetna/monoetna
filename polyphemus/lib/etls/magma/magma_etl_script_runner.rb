require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'
require_relative '../../data_processing/flow_jo_dsl'

class MagmaEtlScriptRunner < EtlScriptRunner
  include FlowJoDsl

  attr_reader :record, :magma_client, :update_request
  def model_name
    if @script_name_parts.length < 3
      raise "Improperly formatted script name #{@file_path}.  Script name should start with '<project name>+<model name>+'."
    end

    @script_name_parts[1]
  end

  def run(record, magma_client:, commit: true)
    # Make the record available as an instance variable named based on the model_name.
    instance_variable_set(:"@#{model_name}", record)
    @record = record
    @magma_client = magma_client
    @update_request = Etna::Clients::Magma::UpdateRequest.new

    run_script

    if commit
      magma_client.update(update_request)
    end
  end
end
