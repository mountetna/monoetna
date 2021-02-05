require 'json'
require_relative './vulcan_controller'

class DataController < Vulcan::Controller
  def workflow_name
    raise Etna::NotFound, "No data for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]
    @params[:workflow_name]
  end

  def cell_hash
    raise Etna::BadRequest, "Relative paths not allowed in data urls" if @params[:cell_hash].include?(".")
    @params[:cell_hash]
  end

  def project_name
    raise Etna::BadRequest, "Relative paths not allowed in data urls" if @params[:project_name].include?(".")
    @params[:project_name]
  end

  def data_filename
    @params[:data_filename]
  end

  def fetch
    return [
        200,
        { 'X-Sendfile' => Vulcan::Storage.data_path(project_name: project_name, cell_hash: cell_hash, data_filename: data_filename),
          'Content-Disposition' => "attachment; filename=#{data_filename}"
        },
        [ '' ]
    ]
  end
end

