require 'etna'

class VulcanV2Controller < Vulcan::Controller

  def workspace_create
    success_json({'it works!': true})
  end

  def pipeline_run
    success_json({'it works!': true})
  end

  def params
    success_json({'it works!': true})
  end

  def list_pipelines
    success_json({'it works!': true})
  end

  def list_workspaces
    success_json({'it works!': true})
  end


end
