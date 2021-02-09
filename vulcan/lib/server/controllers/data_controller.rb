require 'json'
require_relative './vulcan_controller'

class DataController < Vulcan::Controller
  def fetch
    # This is a stub API that will be used for development only.
    # We expect that the CWL YAML file will eventually be served
    #   by Archimedes, not Vulcan itself.

    raise Etna::NotFound, "No data for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    raise Etna::NotFound, "No data for value #{@params[:data]}." unless ["steps", "pools", "umap_data"].include?(@params[:data])

    raise Etna::BadRequest, "Invalid format parameter: #{params[:format]}." if @params[:format] && !["json", "yaml"].include?(@params[:format])

    case @params[:data]
    when "steps"
      if @params[:format] && "json" == @params[:format]
        filename = "steps.json"
        mimetype = "application/json"
      else
        filename = "steps.yaml"
        mimetype = "text/yaml"
      end
    when "pools"
      filename = "pools.json"
      mimetype = "application/json"
    when "umap_data"
      filename = "umap_data.json"
      mimetype = "application/json"
    end
    success(File.read(File.join(
      File.dirname(__FILE__),
      "../data/#{filename}")), mimetype)

  end
end

