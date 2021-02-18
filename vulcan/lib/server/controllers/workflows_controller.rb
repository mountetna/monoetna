require 'json'
require_relative './vulcan_controller'

class WorkflowsController < Vulcan::Controller
  # TODO: When we add real workflow persistent, obviously which workflows are available might change by project.
  # This shares the hard coded projects across all projects from the start.
  def fetch
    success_json({
        'workflows': Dir.glob(Vulcan.instance.config(:workflows_folder) + "/*.{cwl,yaml,yml}").map do |path|
          cwl_name = ::File.basename(path)
          begin
            workflow = Etna::Cwl::Workflow.from_yaml_file(cwl_name)
          rescue => e
            Vulcan.instance.logger.log_error(e)
            next nil
          end

          workflow.as_steps_json(cwl_name)
        end.select { |v| v }
    })
  end

  def submit
    # This is a stub API that will be used for development only.
    # Submit input data and get the specified output status back.

    raise Etna::NotFound, "No data for workflow #{@params[:workflow_name]}." unless "umap" == @params[:workflow_name]

    success(File.read(File.join(
      File.dirname(__FILE__),
      "../data/#{@params[:status]}.json")), 'application/json')

  end
end

