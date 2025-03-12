class Polyphemus
  class Log < Sequel::Model
    def to_report
      {
        application: application,
        project_name: project_name,
        user: user,
        event: event,
        message: message,
        created_at: created_at.iso8601
      }
    end
  end
end
