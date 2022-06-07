class Vulcan
  class WorkflowSnapshot < Sequel::Model
    plugin :timestamps, update_on_create: true

    def self.metadata_params
      [
        :authors,
        :projects,
        :vignette,
        :display_name,
        :description,
        :query_action,
        :input_query_map,
        :tags,
        :image,
      ]
    end

    def self.from_workflow_name(workflow_name:, figure_id:)
      cwl_yaml = Etna::Cwl::Workflow.raw_yaml_from_file(workflow_name)
      metadata = Etna::Cwl::Workflow.metadata(workflow_name)
      Vulcan::WorkflowSnapshot.create({
        figure_id: figure_id,
        cwl_yaml: YAML.dump(cwl_yaml),
      }.update(
        Vulcan::WorkflowSnapshot.snake_case_keys(metadata.slice(*Vulcan::WorkflowSnapshot.metadata_params))
      ))
    end

    def self.from_snapshot(previous_snapshot:, figure_id:)
      Vulcan::WorkflowSnapshot.create({
        figure_id: figure_id,
      }.update(previous_snapshot.to_hash))
    end

    def as_steps_json
      Etna::Cwl::Workflow.from_yaml(cwl_as_yaml).as_steps_json(figure.workflow_name)
    end

    def to_hash
      {}.tap do |result|
        [:cwl_yaml].concat(Vulcan::WorkflowSnapshot.metadata_params).each do |param|
          result[param] = self[param]
        end
      end
    end

    def cwl_as_yaml
      YAML.safe_load(cwl_yaml)
    end

    private

    def figure
      @figure ||= Vulcan::Figure.where(id: figure_id).first
    end

    def self.snake_case_keys(yaml_workflow)
      {}.tap do |snake_case|
        yaml_workflow.each do |key, value|
          snake_case[key.to_s.snake_case] = value
        end
      end
    end
  end
end
