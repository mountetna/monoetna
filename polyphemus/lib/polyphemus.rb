require_relative 'commands'
require 'etna/clients/magma'

class Polyphemus
  include Etna::Application

  def magma_client(token = ENV['TOKEN'])
    Etna::Clients::Magma.new(token: token, host: config(:magma)[:host])
  end
end
