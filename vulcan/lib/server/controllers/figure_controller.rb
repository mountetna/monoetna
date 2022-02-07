require 'json'
require_relative './vulcan_controller'

class FigureController < Vulcan::Controller
  def fetch
    success_json(
      figures: Vulcan::Figure.where(project_name: @params[:project_name]).all.map(&:to_hash)
    )
  end

  def get
    figure = Vulcan::Figure.where(
      project_name: @params[:project_name],
      id: @params[:figure_id]
    ).first

    raise Etna::FileNotFound unless figure

    success_json(figure.to_hash)
  end

  def create
    now = DateTime.now
    figure = Vulcan::Figure.create(
      {
        id: Vulcan::Figure.next_id,
        created_at: now,
        updated_at: now,
      }.update(
        @params.slice(
          :project_name,
          :workflow_name,
          :inputs,
          :title
        )
      )
    )

    success_json(figure.to_hash)
  end
end

