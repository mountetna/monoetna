require_relative 'controller'

class ConfigurationController < Polyphemus::Controller
  def action
    app = Polyphemus.instance
    return_data = {
        app.environment => {
          magma: app.config(:magma)&.dig(:host),
          metis: app.config(:metis)&.dig(:host),
          janus: app.config(:janus)&.dig(:host),
          timur: app.config(:timur)&.dig(:host),
          polyphemus: app.config(:polyphemus)&.dig(:host),
          auth_redirect: app.config(:auth_redirect),
        }
    }
    return success(return_data.to_json, 'application/json')
  end
end
