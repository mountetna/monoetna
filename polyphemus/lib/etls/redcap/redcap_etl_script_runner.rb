require_relative 'lib/client'
require_relative 'lib/model'
require_relative 'lib/loader'
require_relative 'lib/script'
require_relative 'lib/value'
require_relative 'lib/form'

require_relative '../etl_script_runner'
require_relative '../../magma_record_etl'

class Polyphemus
  class RedcapEtlScriptRunner < EtlScriptRunner

    attr_reader :magma_client, :update_request, :model_names, :redcap_token, :redcap_host, :magma_host, :dateshift_salt

    def initialize(project_name:, model_names:, redcap_token:, redcap_host:, magma_host:, dateshift_salt:)
      # Override initialize, user won't be passing in a filename directly.
      @file_path = File.join(File.dirname(__FILE__), 'projects', "#{project_name}.rb")
      @project_name = project_name
      @model_names = model_names
      @redcap_token = redcap_token

      raise "REDCap host must use https://" if redcap_host.start_with?("http://")
      raise "Magma host must use https://" if magma_host.start_with?("http://")

      @redcap_host = redcap_host
      @magma_host = magma_host
      @dateshift_salt = dateshift_salt
    end

    def run(magma_client:, commit: false)
      run_script(self.__binding__)

      loader = Redcap::Loader.new(full_config, magma_client)

      records = loader.run
      puts records

      @update_request = Etna::Clients::Magma::UpdateRequest.new(
        project_name: @project_name,
        revisions: records)

      if commit
        puts "Posting revisions"
        magma_client.update_json(update_request)
      end
      return records
    end

    protected

    def disable_model?(model_name)
      ['all'] == model_names ? false : !model_names.include?(model_name)
    end

    def full_config
      config.update({
        tokens: [ redcap_token ],
        dateshift_salt: dateshift_salt || Polyphemus.instance.sign.uid,
        redcap_host: redcap_host,
        magma_host: magma_host,
        project_name: @project_name,
      })
    end

    def define_model(model_name)
      Object.const_set(model_name, Class.new(Redcap::Model) {})
    end
  end
end
