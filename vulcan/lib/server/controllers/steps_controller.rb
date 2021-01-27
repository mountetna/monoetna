require 'json'
require_relative './vulcan_controller'

class StepsController < Vulcan::Controller
  def fetch
    # This is a stub API that will be used for development only.
    # We expect that the CWL YAML file will eventually be served
    #   by Archimedes, not Vulcan itself.

    raise Etna::NotFound, "No steps for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    success(File.read(File.join(
      File.dirname(__FILE__),
      '../steps/umap.yaml')), 'text/yaml')

  end
end

