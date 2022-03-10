class Vulcan
  class Figure < Sequel::Model
    plugin :timestamps, update_on_create: true

    def self.next_id
      Vulcan.instance.db.get{nextval('figures_ids')}
    end

    def to_hash
      {
        figure_id: figure_id,
        project_name: project_name,
        workflow_name: workflow_name,
        inputs: inputs,
        author: author,
        title: title,
        tags: tags,
        created_at: created_at.iso8601,
        updated_at: updated_at.iso8601
      }
    end
  end
end
