require 'byebug'

class StatsController < Janus::Controller
  def project_stats
    projects_to_include = @params[:projects]

    if projects_to_include.nil?
      projects = Project.all
      permissions = Permission.all
    else
      projects = Project
        .where(project_name: projects_to_include)
        .all
      permissions = Permission
        .left_join(:projects, id: :project_id)
        .where(project_name: projects_to_include)
        .all
    end

    projects = projects.map do |proj|
      {
        project_name: proj.project_name,
        project_name_full: proj.project_name_full,
        resource: proj.resource,
        requires_agreement: proj.requires_agreement,
        principal_investigators: [],
        user_count: 0,
      }
    end

    permissions.each do |perm|
      project_name = perm.project.project_name
      proj = projects.find { |proj| proj[:project_name] == project_name }

      proj[:user_count] += 1

      if perm[:affiliation] == 'PI'

        user = User[perm[:user_id]]

        proj[:principal_investigators].push({
          name: user[:name],
          email: user[:email],
        })
      end
    end

    user_count = User.count

    success_json({projects: projects, user_count: user_count})
  end
end
