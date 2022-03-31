class Vulcan
  class Figure < Sequel::Model
    plugin :timestamps, update_on_create: true

    def self.next_id
      Vulcan.instance.db.get { nextval("figures_ids") }
    end

    def session
      @session ||= Session.from_figure(self)
    end

    def thumbnails(storage:)
      return [] unless storage

      cache = {}

      begin
        workflow = session.workflow
      rescue Etna::Cwl::WorkflowNotFound
        return []
      end

      workflow.outputs.select do |output|
        output.format == "image/png"
      end.map do |output|
        step_name = output.outputSource.first
        output_name = output.outputSource.last

        step_target = session.orchestration.build_target_for(step_name, cache)

        next unless step_target.build_outputs.key?(output_name)

        thumbnail = step_target.build_outputs[output_name]

        next unless ::File.exists?(thumbnail.data_path(storage))

        next storage.data_url(
               project_name: thumbnail.project_name,
               cell_hash: thumbnail.cell_hash,
               data_filename: thumbnail.data_filename,
             )
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
        updated_at: updated_at.iso8601,
      }
    end
  end
end
