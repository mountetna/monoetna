require 'json'
require_relative './vulcan_controller'

class FigureController < Vulcan::Controller
  def fetch
    success_json(
      figures: Vulcan::Figure.where(
        project_name: @params[:project_name],
        archived: false).all.map do |f|
        f.to_hash(storage: storage)
      end
    )
  end

  def get
    figure = Vulcan::Figure[
      @params[:project_name],
      @params[:figure_id]
    ]

    raise Etna::NotFound unless figure

    success_json(figure.to_hash(storage: storage))
  end

  def create
    now = DateTime.now
    figure = Vulcan::Figure.from_payload(
      {
        figure_id: Vulcan::Figure.next_id,
        author: @user.name,
        created_at: now,
        updated_at: now,
      }.update(
        @params.slice(*figure_params)
      ).update(
        dependencies: dependency_shas
      )
    )
    success_json(figure.to_hash)
  rescue ArgumentError => e
    failure(422, e.message)
  end

  def update
    require_params(:comment)

    figure = Vulcan::Figure[
      @params[:project_name],
      @params[:figure_id]
    ]

    raise Etna::NotFound unless figure

    now = DateTime.now

    new_figure_data = {
      figure_id: figure.figure_id,
      author: @user.name,
      created_at: figure.created_at,
      updated_at: now,
      comment: @params[:comment]
    }.update(
      figure.to_hash.slice(*figure_params)
    ).update(
      @params.slice(*figure_params)
    )

    new_figure = Vulcan::Figure.create(new_figure_data)

    # Only archive the original figure once the new figure
    #   exists, otherwise we risk losing the figure if
    #   an exception gets thrown during figure creation.
    figure.modified!(:updated_at)
    figure.update(archived: true)
  
    begin
      new_figure.update_dependencies
      new_figure.take_snapshot
    end if @params[:update_dependencies]

    success_json(new_figure.to_hash)
  rescue ArgumentError => e
    failure(422, e.message)
  end

  def revisions
    figures = Vulcan::Figure.where(
      project_name: @params[:project_name],
      figure_id: @params[:figure_id]
    ).reverse_order(:updated_at).all.sort_by{|e| e.archived ? 1 : 0 }

    success_json(figures.map(&:to_revision))
  end

  def delete
    figure = Vulcan::Figure[
      @params[:project_name],
      @params[:figure_id]
    ]

    raise Etna::NotFound unless figure

    figure.update(
      archived: true
    )

    success_json(figure.to_hash)
  end

  private

  def figure_params
    [
      :project_name,
      :workflow_name,
      :inputs,
      :title,
      :tags,
      :dependencies
    ]
  end

  def dependency_shas
    Vulcan.instance.dependency_manager.dependency_shas.to_json
  end
end

