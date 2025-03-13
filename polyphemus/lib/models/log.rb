class Polyphemus
  class Log < Sequel::Model
    def to_report
      {
        id: id,
        application: application,
        project_name: project_name,
        user: user,
        event: event,
        message: message,
        payload: !!payload,
        created_at: created_at.iso8601
      }
    end
  end
end
