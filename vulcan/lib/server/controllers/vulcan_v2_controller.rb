require 'etna'

class VulcanV2Controller < Vulcan::Controller

  def init
    success_json({'it works!': true})
  end

  def list
    success_json({'it works!': true})
  end
end
