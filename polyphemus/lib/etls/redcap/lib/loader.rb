module Redcap
  class Loader
    def initialize(config, magma_client)
      @config = config
      @magma_client = magma_client
    end

    def run
      magma_models = @magma_client.retrieve(
        Etna::Clients::Magma::RetrievalRequest.new(project_name: @config[:project_name])
      ).models

      records = {}
      models.each do |model_name, model_config|
        next if model_config[:disabled]
        records[ model_name ] = {}

        tokens.each do |token|
          client = Redcap::Client.new(@config[:redcap_host], token)

          records[ model_name ].update(
            Redcap::Model.create(
              model_name,
              model_config[:scripts],
              magma_models.model(model_name.to_s).template,
              dateshift_salt
            ).load(client)
          )
        end
      end

      patch_tables(records, magma_models)

      records
    end

    private

    def patch_tables(records, magma_models)
      # find any table attributes
      magma_models.model_keys.each do |model_name|
        magma_models.model(model_name.to_s).template.attributes.tap do |atts|
          atts.attribute_keys.each do |att_name|
            if atts.attribute(att_name).attribute_type == 'table'
              # look for corresponding entries in the database
              link_model_name = atts.attribute(att_name).link_model_name.to_sym

              if records[link_model_name]
                records[link_model_name].to_a.group_by do |(record_name,record)|
                  record[ model_name.to_sym ]
                end.each do |model_record_name, revisions|
                  temp_id_names = revisions.map(&:first)
                  records[ model_name.to_sym ] ||= {}
                  records[ model_name.to_sym ][ model_record_name ] ||= {}
                  records[ model_name.to_sym ][ model_record_name ][ att_name ] = temp_id_names
                end
              end
            end
          end
        end
      end
    end

    def models
      @config[:models]
    end

    def dateshift_salt
      @config[:dateshift_salt]
    end

    def tokens
      @config[:tokens]
    end
  end
end
