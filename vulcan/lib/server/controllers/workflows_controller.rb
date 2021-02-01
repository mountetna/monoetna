require 'json'
require_relative './vulcan_controller'

class WorkflowsController < Vulcan::Controller
  def fetch
    success_json({
      umap: JSON.load(
          File.read(File.join(
            File.dirname(__FILE__),
            '../workflows/umap.json')),
        symbolize_names: true)
    })
  end
end

