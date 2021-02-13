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

  def data_path
    @data_path ||= begin
      storage.data_path(project_name: project_name, cell_hash: cell_hash, data_filename: data_filename).tap do |path|
        unless ::File.exists?(path)
          raise Etna::NotFound, "File #{data_filename} for hash #{cell_hash} in project #{project_name} was not found"
        end
      end
    end
  end

  # TODO: Add content type and additional metadata from a cell's payload.
  def fetch
    return [
        200,
        { 'X-Sendfile' => data_path,
          'Content-Disposition' => "attachment; filename=#{data_filename}"
        },
        [ '' ]
    ]
  end
end

