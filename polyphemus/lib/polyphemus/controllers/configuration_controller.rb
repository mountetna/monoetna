require_relative 'controller'

class ConfigurationController < Polyphemus::Controller
  def action
    return failure(401, errors: [ 'You are unauthorized' ]) unless @user && @user.can_view?(@project_name)

    return_data = {
      magma: Polyphemus.instance.config(:magma)[:host],
      metis: Polyphemus.instance.config(:metis)[:host],
      janus: Polyphemus.instance.config(:janus)[:host],
      timur: Polyphemus.instance.config(:timur)[:host],
      polyphemus: Polyphemus.instance.config(:polyphemus)[:host],
    }
    return success(return_data.to_json, 'application/json')
  end
end
