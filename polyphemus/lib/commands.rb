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

  class GetTree < Etna::Command
    usage 'Open a console with a connected Polyphemus instance.'

    def execute
      client = Polyphemus.instance.magma_client
      models = client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: 'mvir1')).models
      p models.to_directed_graph(true).descendants('patient')
    end

    def setup(config)
      super
    end
  end

  class GetMetisFolders < Etna::Command
    usage "Fetch a list of Metis folders from a bucket"

    def execute
      client = Polyphemus.instance.metis_client
      folders = client.folder_list_all_folders(
        Etna::Clients::Metis::ListFoldersRequest.new(project_name: 'mvir1', bucket_name: 'data')).folders
      folders.all.each { |f| p f.folder_path }
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
