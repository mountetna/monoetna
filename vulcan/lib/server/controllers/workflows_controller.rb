require 'json'
require_relative './vulcan_controller'

class WorkflowsController < Vulcan::Controller
  # TODO: When we add real workflow persistent, obviously which workflows are available might change by project.
  # This shares the hard coded projects across all projects from the start.
  def fetch
    success_json({
        'workflows': Dir.glob(Vulcan.instance.config(:workflows_folder) + "/*.{cwl,yaml,yml}").sort_by { |f| ::File.basename(f) }.map do |path|
          cwl_name = ::File.basename(path)
          begin
            workflow = Etna::Cwl::Workflow.from_yaml_file(cwl_name)
          rescue => e
            Vulcan.instance.logger.log_error(e)
            next nil
          end

          # Really we should fetch the metadata from a more
          #   permanent store. Right now storing them in
          #   JSON files in the workflows folder.
          metadata_file = Vulcan.instance.config(:workflows_folder) + "/#{cwl_name.gsub('.cwl', '')}.metadata.json"

          workflow.as_steps_json(cwl_name).update(
            File.exists?(metadata_file) ?
              JSON.parse(File.read(metadata_file),
              symbolize_names: true) :
              {}).tap do |workflow_json|
            # Workflows without specific project restrictions are assumed to be available to all the user's projects.
            workflow_json[:projects] ||= @user.projects.map(&:to_s)
            workflow_json[:projects] = workflow_json[:projects].select { |p| @user.can_view?(p) }
          end
        end.compact.select do |v|
          # We only want workflows where the user is
          #   authorized for at least one listed projects
          !v[:projects].empty?
        end
    })
  end
rescue => e
  Vulcan.instance.logger.log_error(e)
  raise Etna::BadRequest.new(e.message)
end

