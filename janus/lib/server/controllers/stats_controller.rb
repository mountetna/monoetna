class StatsController < Janus::Controller
  def stats
    projects = Project.all.map do |proj|
      {
        project_name: proj.project_name,
        project_name_full: proj.project_name_full,
        resource: proj.resource,
      }
    end

    Permission.all.each do |perm|
      project_name = perm.project.project_name

      projects[project_name][:user_count] = 0 unless projects[project_name].key?(:user_count)
      projects[project_name][:user_count] += 1
    end

    user_count = User.all.count

    success_json({projects: projects, user_count: user_count})
  end
end
