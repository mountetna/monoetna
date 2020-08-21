require_relative 'commands'
require 'etna/clients/magma'

class Polyphemus
  include Etna::Application

  def magma_client
    @magma_client ||= Etna::Clients::Magma.new({ token: ENV['TOKEN'] }.update(config(:magma)))
  end
end
