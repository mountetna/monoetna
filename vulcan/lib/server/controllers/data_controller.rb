require 'json'
require_relative './vulcan_controller'

class DataController < Vulcan::Controller
  def fetch
    # This is a stub API that will be used for development only.
    # We expect that the CWL YAML file will eventually be served
    #   by Archimedes, not Vulcan itself.

    raise Etna::NotFound, "No data for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    raise Etna::NotFound, "No data for value #{@params[:data]}." unless ["steps", "pools"].include?(@params[:data])

    case @params[:data]
    when "steps"
      filename = "steps.yaml"
      mimetype = "text/yaml"
    when "pools"
      filename = "pools.json"
      mimetype = "application/json"
    end
    success(File.read(File.join(
      File.dirname(__FILE__),
      "../data/#{filename}")), mimetype)

  end
end

