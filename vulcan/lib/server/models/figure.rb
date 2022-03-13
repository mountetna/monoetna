class Vulcan
  class Figure < Sequel::Model
    plugin :timestamps, update_on_create: true

    def self.next_id
      Vulcan.instance.db.get{nextval('figures_ids')}
    end

    def session
      @session ||= Session.from_figure(self)
    end

    def thumbnails(storage:)
      return [] unless storage
      
      cache = {}

      session.workflow.steps.map do |step|
        case step.ui_output_name
        when 'plotly.cwl'
          step_target = session.orchestration.build_target_for(step.id, cache)

          next unless step_target.input_files.length > 1

          thumbnail = step_target.input_files[1]

          next unless ::File.exists?(thumbnail.data_path(storage))

          next storage.data_url(
            project_name: thumbnail.project_name,
            cell_hash: thumbnail.cell_hash,
            data_filename: thumbnail.data_filename
          )
        else
          next
        end
      end.compact
    end

    def to_hash(storage: nil)
      {
        figure_id: figure_id,
        project_name: project_name,
        workflow_name: workflow_name,
        key: "#{project_name}/#{figure_id}",
        inputs: inputs.to_hash,
        author: author,
        title: title,
        tags: tags,
        thumbnails: thumbnails(storage: storage),
        created_at: created_at.iso8601,
        updated_at: updated_at.iso8601
      }
    end
  end
end
