# require_relative 'controller'

# class ConfigurationController < Polyphemus::Controller
#   def action
#     app = Polyphemus.instance
#     return_data = {
#         app.environment => {
#           docker: app.config(:docker)&.slice(:default_tag),
#           magma: app.config(:magma)&.slice(:host),
#           metis: app.config(:metis)&.slice(:host),
#           janus: app.config(:janus)&.slice(:host),
#           timur: app.config(:timur)&.slice(:host),
#           polyphemus: app.config(:polyphemus)&.slice(:host),
#           auth_redirect: app.config(:auth_redirect),
#         }
#     }
#     return success(return_data.to_json, 'application/json')
#   end
# end
