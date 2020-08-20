require_relative 'commands'
require 'etna/clients/magma'
require 'etna/clients/metis'

class Polyphemus
  include Etna::Application

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new({ token: ENV['TOKEN'] }.update(config(:magma)))
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new({ token: ENV['TOKEN'] }.update(config(:metis)))
  end
end
