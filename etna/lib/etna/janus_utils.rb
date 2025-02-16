# Utility class to work with Janus for authz and authn
module Etna
  class JanusUtils
    def initialize
    end

    def projects(token)
      return [] unless has_janus_config?

      janus_client(token).get_projects.projects
    rescue
      # If encounter any issue with Janus, we'll return no projects
      []
    end

    def resource_projects(token)
      projects(token).select do |project|
        project.project_type == 'resource'
      end
    end

    def community_projects(token)
      projects(token).select do |project|
        project_type == 'community'
      end
    end

    def janus_client(token)
      Etna::Clients::Janus.new(
        token: token,
        host: application.config(:janus)[:host],
      )
    end

    def valid_task_token?(token)
      return false unless has_janus_config?

      response = janus_client(token).validate_task_token

      return false unless response.code == "200"

      return true
    end

    private

    def application
      @application ||= Etna::Application.instance
    end

    def has_janus_config?
      application.config(:janus) && application.config(:janus)[:host]
    end
  end
end
