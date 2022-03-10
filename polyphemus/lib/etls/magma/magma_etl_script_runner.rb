require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'
require_relative '../../data_processing/flow_jo_dsl'

class Polyphemus
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
      @update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project_name)

      # In production, seems like we can't pass
      #   self.__binding__ here -- throws an UndefinedMethod exception.
      run_script(self.get_binding)

      if commit
        magma_client.update_json(update_request)
      end
    end

    def get_binding
      binding
    end
  end

  # TODO: Nuke all of this, fold the ipi flojo processing into airflow instead.
  class MagmaRecordScriptEtl < Polyphemus::MagmaRecordEtl
    def initialize(**args)
      super(**{project_model_pairs: [[script_runner.project_name, script_runner.model_name]]}.update(args))
    end

    def script_runner
      self.class.instance_variable_get(:@script_runner)
    end

    def process(cursor, batch)
      batch.each do |record|
        script_runner.run(record, magma_client: magma_client)
      end
    end

    module Scripts
      # For each script file, do the thing.
      Dir[::File.join(__dir__, "scripts/*.rb")].each do |file|
        script_runner = MagmaEtlScriptRunner.new(file)
        name = script_runner.script_name_parts.map(&:capitalize).join('')
        const_set(:"#{name}", Class.new(MagmaRecordScriptEtl) do
          instance_variable_set(:@script_runner, script_runner)
        end)
      end
    end
  end
end