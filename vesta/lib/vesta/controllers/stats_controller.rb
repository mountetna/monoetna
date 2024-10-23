require_relative 'controller'

class StatsController < Vesta::Controller
  def global_stats
    # TODO: include filter by time

    stats = Vesta::GlobalStats.all

    success_json(stats.map(&:to_hash))
  end

  def project_stats
    projects_to_include = @params[:projects]

    stats = Vesta::ProjectStats
    unless projects_to_include.nil?
      stats = stats.where(name: projects_to_include)
    end

    success_json(stats.map(&:to_hash))
  end
end
