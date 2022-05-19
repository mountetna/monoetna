class Vulcan
  class Figure < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_one :workflow_snapshot

    def self.from_reference_id(reference_id)
      Vulcan::Figure.where(id: reference_id).first
    end

    def self.next_id
      Vulcan.instance.db.get { nextval("figures_ids") }
    end

    def self.[](project_name, figure_id)
      Vulcan::Figure.where(
        project_name: project_name,
        figure_id: figure_id,
        archived: false
      ).first
    end

    def after_save
      # Use this hook to also save the Vulcan workflow snapshot
      # Only when the workflow has changed and no current snapshot for the figure
      begin
        previous_figure = Vulcan::Figure.where(figure_id: self.figure_id).exclude(id: self.id).last
        previous_snapshot = Vulcan::WorkflowSnapshot.where(figure_id: previous_figure.id).first if previous_figure
      
        current_workflow = current_workflow_json
      
        Vulcan::WorkflowSnapshot.from_workflow_json(
          figure_id: self.id,
          json_workflow: workflow_changed?(previous_snapshot, current_workflow) && previous_snapshot ?
            previous_snapshot.to_workflow_json :
            current_workflow
        )
    
        self.refresh
      end if !has_snapshot?
      
      super
    end

    def has_snapshot?
      !!Vulcan::WorkflowSnapshot.where(figure_id: self.id).first
    end

    def workflow_changed?(previous_snapshot, current_workflow)
      previous_snapshot&.to_workflow_json != current_workflow
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

    def to_revision
      {
        inputs: inputs.to_hash,
        title: title,
        tags: tags,
        updated_at: updated_at.iso8601,
        comment: comment
      }
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
        dependencies: dependencies,
        workflow_snapshot: workflow_snapshot&.to_workflow_json
      }
    end

    private

    def current_workflow_json
      workflow = session.workflow
      json_workflow = workflow.as_steps_json(self.workflow_name)
      json_workflow.update(Etna::Cwl::Workflow.metadata(json_workflow[:name]))
    end
  end
end
