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
      figure_id: @params[:figure_id]
    ).order(:updated_at).first

    raise Etna::NotFound unless figure

    success_json(figure.to_hash)
  end

  def create
    now = DateTime.now
    figure = Vulcan::Figure.create(
      {
        figure_id: Vulcan::Figure.next_id,
        author: @user.name,
        created_at: now,
        updated_at: now,
      }.update(
        @params.slice(*figure_params)
      )
    )
    success_json(figure.to_hash)
  end

  def update
    figure = Vulcan::Figure.where(
      project_name: @params[:project_name],
      figure_id: @params[:figure_id]
    ).first

    raise Etna::NotFound unless figure

    figure.update(
        @params.slice(*figure_params)
    )

    success_json(figure.to_hash)
  end

  def delete
    figure = Vulcan::Figure.where(
      project_name: @params[:project_name],
      figure_id: @params[:figure_id]
    ).first

    raise Etna::NotFound unless figure

    figure.delete

    success_json(figure.to_hash)
  end

  private

  def figure_params
    [
      :project_name,
      :workflow_name,
      :inputs,
      :title,
      :tags
    ]
  end
end

