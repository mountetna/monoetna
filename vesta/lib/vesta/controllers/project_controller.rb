require_relative 'controller'

class ProjectController < Vesta::Controller
  def projects
    projects_to_include = @params[:projects]

    projects = Vesta::Project
    projects = projects.where(name: projects_to_include) unless projects_to_include.nil?

    success_json(projects.map(&:to_hash))
  end
end
