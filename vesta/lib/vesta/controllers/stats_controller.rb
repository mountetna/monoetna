require_relative 'controller'

class StatsController < Vesta::Controller
  def global_stats
    # TODO: include filter by time

    stats = Vesta::GlobalStats.all

    success_json(stats)
  end

  def project_stats
    projects_to_include = @params[:projects]

    stats = Vesta::ProjectStats
    stats = stats.where(name: projects_to_include) unless projects_to_include.nil?

    success_json(stats.all)
  end
end
