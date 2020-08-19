require 'date'
require 'logger'

class Polyphemus
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Polyphemus.instance.commands.each do |name,cmd|
        puts cmd.usage
      end
    end
  end

  class CascadeMvirPatientWaiverToRestricted < Etna::Command
    usage 'Updates any models whose restricted access does not match its relationship to patient'

    def execute
      client = Polyphemus.instance.magma_client(Polyphemus.instance.config(:mvir1)[:token])
      models = client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: 'mvir1')).models

      client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(model_name: 'patient', attribute_names: ['identifier', '']))

      models.to_directed_graph(true).descendants('patient').keys.dup.push('patient').each do |model_name|
        model_name
      end
    end
  end

  class AddMvirRestrictedColumns < Etna::Command
    usage 'Adds restricted columns to all models descending from patient that currently do not have it'

    def execute
      client = Polyphemus.instance.magma_client(Polyphemus.instance.config(:mvir1)[:token])

      models = client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: 'mvir1')).models
      request = Etna::Clients::Magma::UpdateModelRequest.new(project_name: 'mvir1')
      models.to_directed_graph(true).descendants('patient').keys.dup.push('patient').each do |model_name|
        next if models.model(model_name).attributes.attribute_keys.include?('restricted')
        request.add_action(Etna::Clients::Magma::AddAttributeAction.new(
            model_name: model_name,
            attribute_name: 'restricted',
            type: Etna::Clients::Magma::AttributeType::BOOLEAN.to_s,
            restricted: true,
            description: 'Controls access to the data, mostly as a means of COMET waiver consent status.'
        ))
      end
      client.update_model(request)
    end

    def setup(config)
      super
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Polyphemus instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
    end
  end
end
