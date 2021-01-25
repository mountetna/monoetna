require 'json'

class StepsController < Vulcan::Controller
  def fetch
    raise Etna::NotFound, "No steps for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    steps = JSON.parse(
      File.expand_path('../steps/umap.json', __FILE__),
      symbolize_names: true)

    success_json(steps)
  end
end

