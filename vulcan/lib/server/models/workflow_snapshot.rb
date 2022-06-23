require "active_support/inflector"

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
        Vulcan::WorkflowSnapshot.snake_case_keys(metadata).slice(*Vulcan::WorkflowSnapshot.metadata_params)
      ).update(
        scripts: Etna::Cwl::Workflow.step_scripts(cwl_yaml),
      ))
    end

    def self.from_snapshot(previous_snapshot:, figure_id:)
      Vulcan::WorkflowSnapshot.create({
        figure_id: figure_id,
      }.update(previous_snapshot.to_hash))
    end

    def as_steps_json_w_metadata
      Etna::Cwl::Workflow.from_yaml(cwl_as_yaml).as_steps_json(figure.workflow_name).update(
        params_to_hash(Vulcan::WorkflowSnapshot.metadata_params, camel_case: true)
      )
    end

    def cwl_as_yaml
      YAML.safe_load(cwl_yaml)
    end

    def to_hash
      params_to_hash([:cwl_yaml, :scripts].concat(Vulcan::WorkflowSnapshot.metadata_params))
    end

    private

    def params_to_hash(params, camel_case: false)
      {}.tap do |result|
        params.each do |param|
          if camel_case
            result[param.to_s.camelize(:lower).to_sym] = self[param]
          else
            result[param] = self[param]
          end
        end
      end
    end

    def figure
      @figure ||= Vulcan::Figure.where(id: figure_id).first
    end

    def self.snake_case_keys(camel_case_hash)
      {}.tap do |snake_case|
        camel_case_hash.each do |key, value|
          snake_case[key.to_s.snake_case.to_sym] = value
        end
      end
    end
  end
end
