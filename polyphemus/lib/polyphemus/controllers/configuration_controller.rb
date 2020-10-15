require_relative 'controller'

class ConfigurationController < Polyphemus::Controller
  def action
    puts "HERE I AM"
    return_data = {
      magma: Polyphemus.instance.config(:magma)&.dig(:host),
      metis: Polyphemus.instance.config(:metis)&.dig(:host),
      janus: Polyphemus.instance.config(:janus)&.dig(:host),
      timur: Polyphemus.instance.config(:timur)&.dig(:host),
      polyphemus: Polyphemus.instance.config(:polyphemus)&.dig(:host)
    }
    return success(return_data.to_json, 'application/json')
  end
end
