require 'json'
require_relative './vulcan_controller'

class FigureController < Vulcan::Controller
  # TODO: When we add real workflow persistent, obviously which workflows are available might change by project.
  # This shares the hard coded projects across all projects from the start.
  def fetch
    success_json(
      figures: Vulcan::Figures.where(project_name: @params[:project_name]).all.map(&:to_hash)
    )
  end
end

