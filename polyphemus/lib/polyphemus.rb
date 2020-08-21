require_relative 'commands'
require 'etna/clients/magma'
require 'etna/clients/metis'

class Polyphemus
  include Etna::Application

  def magma_client(token = ENV['TOKEN'])
    Etna::Clients::Magma.new(token: token, host: config(:magma)[:host])
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new({ token: ENV['TOKEN'] }.update(config(:metis)))
  end
end
