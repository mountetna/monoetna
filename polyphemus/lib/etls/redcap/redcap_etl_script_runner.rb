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

    attr_reader :magma_client, :update_request, :model_names, :redcap_token

    def initialize(project_name:, model_names:, redcap_token:)
      # Override initialize, user won't be passing in a filename directly.
      # Path is relative to where this is invoked...TODO: need to make this work for
      #   the controller path.
      @file_path = "lib/etls/redcap/projects/#{project_name}.rb"
      @project_name = project_name
      @model_names = model_names
      @redcap_token = redcap_token
    end

    def run(magma_client:, commit: false)
      run_script(self.__binding__)

      loader = Redcap::Loader.new(config, magma_client)

      records = loader.run
      puts records

      @update_request = Etna::Clients::Magma::UpdateRequest.new(
        project_name: @project_name,
        revisions: records)

      puts "Posting revisions"

      # if commit
      #   magma_client.update_json(update_request)
      # end
      return records
    end

    protected

    def disable_model?(model_name)
      ['all'] == model_names ? false : !model_names.include?(model_name)
    end
  end
end
