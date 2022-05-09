require 'json'
require 'net/http'
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
    figure = Vulcan::Figure.create(
      {
        figure_id: Vulcan::Figure.next_id,
        author: @user.name,
        created_at: now,
        updated_at: now,
      }.update(
        @params.slice(*figure_params)
      ).update(
        dependencies: dependency_shas.to_json
      )
    )
    success_json(figure.to_hash)
  end

  def update
    require_params(:comment)

    figure = Vulcan::Figure[
      @params[:project_name],
      @params[:figure_id]
    ]

    raise Etna::NotFound unless figure

    figure.modified!(:updated_at)
    figure.update(archived: true)

    now = DateTime.now

    new_figure = Vulcan::Figure.create(
      {
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
    )

    if @params[:update_dependencies]
      new_figure.update(
        dependencies: dependency_shas.to_json
      )
    end

    success_json(new_figure.to_hash)
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

  def image_sha(image_name:, tag: 'production', namespace: 'etnaagent', registry: 'registry.hub.docker.com/v2/repositories')
    uri = URI("https://#{registry}/#{namespace}/#{image_name}/tags/#{tag}")
    metadata = JSON.parse(Net::HTTP.get(uri), symbolize_names: true)
    metadata[:images].first[:digest]
  end

  def dependency_shas
    Vulcan.instance.dependencies.map do |dependency|
      [dependency, image_sha(image_name: dependency)]
    end.to_h
  end
end

