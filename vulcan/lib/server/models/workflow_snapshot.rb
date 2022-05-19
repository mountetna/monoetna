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
        figure_id: figure_id
      }.update(
        Vulcan::WorkflowSnapshot.snake_case_keys(json_workflow.slice(*Vulcan::WorkflowSnapshot.params))
      ).update(
        other_metadata: json_workflow.slice(*(json_workflow.keys - Vulcan::WorkflowSnapshot.params))
      ))
    end

    def to_workflow_json
      JSON.parse({
        class: "Workflow",
      }.update(Vulcan::WorkflowSnapshot.params.map do |param|
        [param, self[param.to_s.snake_case.to_sym]]
      end.to_h).update(self.other_metadata).to_json)
    end

    private
    
    def self.snake_case_keys(json)
      {}.tap do |snake_case|
        json.keys.each do |key|
          snake_case[key.to_s.snake_case] = json[key]
        end
      end
    end
  end
end
