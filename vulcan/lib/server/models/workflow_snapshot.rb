class Vulcan
  class WorkflowSnapshot < Sequel::Model
    plugin :timestamps, update_on_create: true

    def self.params
      [
        :authors,
        :name,
        :projects,
        :vignette,
        :displayName,
        :description,
        :cwlVersion,
        :inputs,
        :outputs,
        :steps,
        :tags,
        :image,
      ]
    end

    def self.from_workflow_json(json_workflow:, figure_id:)
      Vulcan::WorkflowSnapshot.create({
        figure_id: figure_id,
      }.update(
        Vulcan::WorkflowSnapshot.snake_case_keys(json_workflow.slice(*Vulcan::WorkflowSnapshot.params))
      ).update(
        other_metadata: json_workflow.slice(*(json_workflow.keys - Vulcan::WorkflowSnapshot.params)),
      ))
    end

    def to_workflow_json(symbolize_names: true, for_json_loader: false)
      JSON.parse({
        class: "Workflow",
      }.update(Vulcan::WorkflowSnapshot.params.map do |param|
        value = self[param.to_s.snake_case.to_sym]
        if param == :steps
          [param, [value.first&.map do |step|
            for_json_loader ? step : inject_step_name(step)
          end].compact]
        else
          [param, value]
        end
      end.to_h).update(self.other_metadata).to_json, symbolize_names: symbolize_names)
    end

    private

    def self.snake_case_keys(json_workflow)
      {}.tap do |snake_case|
        json_workflow.each do |key, value|
          if key == :steps
            # Steps are nested one level deep
            snake_case[key.to_s.snake_case] = [value.first&.map do |step|
              Vulcan::WorkflowSnapshot.inject_step_id(step)
            end].compact
          else
            snake_case[key.to_s.snake_case] = value
          end
        end
      end
    end

    def self.inject_step_id(step_json)
      step_json.map do |key, value|
        key.to_sym == :name ? [:id, value] : [key, value]
      end.to_h
    end

    def inject_step_name(step_json)
      step_json.map do |key, value|
        key.to_sym == :id ? [:name, value] : [key, value]
      end.to_h
    end
  end
end
