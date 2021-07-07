require 'json'

require_relative 'lib/client'
require_relative 'lib/model'
require_relative 'lib/template'
require_relative 'lib/loader'
require_relative 'lib/script'
require_relative 'lib/entity'
require_relative 'lib/value'
require_relative 'lib/project'
require_relative 'lib/magma_models'

require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'

class Polyphemus
  class RedcapEtlScriptRunner < EtlScriptRunner

    attr_reader :magma_client, :update_request, :model_names, :redcap_tokens, :redcap_host, :magma_host, :dateshift_salt, :mode

    # Override initialize, user won't be passing in a filename directly as with other ETLs.
    def initialize(project_name:, model_names: "all", redcap_tokens:, redcap_host:, magma_host:, dateshift_salt:, mode: nil)
      raise "No dateshift_salt provided, please provide one." unless dateshift_salt
      raise "Mode must be nil, \"existing\", or \"strict\"." unless [nil, "existing", "strict"].include?(mode)
      raise "Must provide at least one REDCap token." unless redcap_tokens&.length > 0

      @file_path = File.join(File.dirname(__FILE__), 'projects', "#{project_name}.rb")

      raise "Project configuration does not exist." unless File.file?(@file_path)

      @project_name = project_name
      @model_names = model_names
      @redcap_tokens = redcap_tokens

      raise "REDCap host must use https://" unless redcap_host.start_with?("https://")
      raise "Magma host must use https://" unless magma_host.start_with?("https://")

      @redcap_host = redcap_host
      @magma_host = magma_host
      @dateshift_salt = dateshift_salt
      @mode = mode # operating mode: nil, "strict", "existing"
    end

    def run(magma_client:, commit: false, logger: STDOUT)
      yield logger if block_given?

      # For some reason in the Puma environment, can't pass
      #   self.__binding__ here -- throws an UndefinedMethod exception.
      run_script(self.get_binding)

      loader = Redcap::Loader.new(config.update(system_config), @project_name, magma_client, logger)

      all_records, records_to_blank = loader.run

      logger.write(JSON.pretty_generate(all_records))
      logger.write("\n")

      if commit
        logger.write("Posting revisions.\n")

        update_request = Etna::Clients::Magma::UpdateRequest.new(
          project_name: @project_name,
          revisions: all_records)
        magma_client.update_json(update_request)

        logger.write("Revisions saved to Magma.\n")
      end

      summarize(
        logger: logger,
        all_records: all_records,
        records_to_blank: records_to_blank,
        commit: commit)

      return all_records
    rescue => e
      logger.write("#{e.message}\n#{e.backtrace}")
      raise
    end

    def system_config
      {
        tokens: redcap_tokens,
        dateshift_salt: dateshift_salt,
        redcap_host: redcap_host,
        magma_host: magma_host,
        project_name: @project_name,
        models_to_build: model_names,
        mode: mode
      }
    end

    def get_binding
      binding
    end

    protected

    def summarize(logger:, all_records:, records_to_blank:, commit:)
      logger.write(<<-EOM
===============================
Summary of upload
===============================
Project: #{@project_name}
Models: #{model_names}
Mode setting: #{mode}
Committed to Magma: #{commit}
EOM
      )
      all_records.keys.each do |model_name|
        logger.write("#{model_name} records updated: #{all_records[model_name].keys.length}\n")
      end

      if records_to_blank
        records_to_blank.keys.each do |model_name|
          logger.write("#{model_name} records blanked: #{records_to_blank[model_name].length}\n")
        end
      end

      logger.write("===============================\n")
    end

    def define_model(model_name, &block)
      return Kernel.const_get(model_name) if Kernel.const_defined?(model_name)

      # Set some default methods for each model
      Kernel.const_set(model_name, Class.new(Redcap::Model) {
        def identifier(*record_id)
          [
              "::temp", *record_id, rand(36**8).to_s(36)
          ].compact.join('-')
        end

        def patch(id, record)
        end
      })
    end
  end
end
