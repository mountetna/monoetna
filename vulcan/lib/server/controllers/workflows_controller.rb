require 'json'

class WorkflowsController < Vulcan::Controller
  def fetch
    workflows = {
      umap: JSON.parse(
        File.expand_path('../workflows/umap.json', __FILE__),
        symbolize_names: true)
    }

    success_json(workflows)
  end
end

