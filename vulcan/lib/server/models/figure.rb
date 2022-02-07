class Vulcan
  class Figure < Sequel::Model
    def self.next_id
      Vulcan.instance.db.get{nextval('figures_ids')}
    end

    def to_hash
      {
        id: id,
        project_name: project_name,
        workflow_name: workflow_name,
        inputs: inputs,
        title: title,
        created_at: created_at.iso8601,
        updated_at: updated_at.iso8601
      }
    end
  end
end
