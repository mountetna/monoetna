require 'json'
require_relative './vulcan_controller'

class WorkflowsController < Vulcan::Controller
  def fetch
    ::File.dirname
  end

  def submit
    # This is a stub API that will be used for development only.
    # Submit input data and get the specified output status back.
    raise Etna::NotFound, "No data for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    success(File.read(File.join(
      File.dirname(__FILE__),
      "../../../data/#{@params[:status]}.json")), 'application/json')

  end
end

