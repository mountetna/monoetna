require 'json'
require_relative './vulcan_controller'

class StepsController < Vulcan::Controller
  def fetch
    raise Etna::NotFound, "No steps for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    success_json({
      steps: JSON.parse(
        File.read(File.join(
          File.dirname(__FILE__),
          '../steps/umap.json')),
        symbolize_names: true)
    })
  end
end

