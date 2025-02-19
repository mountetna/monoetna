require_relative '../helpers'
require_relative './jobs/etl_job'
require_relative './jobs/metis_linker'
require_relative './jobs/redcap_loader'
require_relative './jobs/sftp_file_discovery'
require_relative './jobs/sftp_deposit_uploader'
require_relative './jobs/sftp_metis_uploader'

class Polyphemus
  module ArgoWorkflowManager
    extend WithLogger

    def submit_workflow(config)
      # Build command with proper escaping
      manifest = Polyphemus::WorkflowManifest.from_workflow_name(config.workflow_type)
      workflow_path = manifest.workflow_path.shellescape
      cmd = [
        "argo submit",
        workflow_path,
        "-p config_id=#{config.config_id}",
        "-p version_number=#{config.version_number}", 
        "-n argo",
        "-o yaml",
        "| grep -m2 -E 'name:|uid:'",
        "| sed 's/name:/Workflow Name:/'",
        "| sed 's/uid:/Workflow UID:/'"
      ]

      logger.info("Submitting workflow #{config.workflow_name}, for project: #{config.project_name}, workflow_type: #{config.workflow_type}, config_id: #{config.config_id}...")
      logger.info("Submitting Argo workflow with command: #{cmd.join(" ")}")
      stdout, stderr, status = Open3.capture3(cmd.join(" "))

      if status.success?
        workflow_name = stdout.match(/Workflow Name: (.+)/)&.[](1)
        workflow_uid = stdout.match(/Workflow UID: (.+)/)&.[](1)
        logger.info("Argo workflow submitted successfully - Name: #{workflow_name}, UID: #{workflow_uid}")
        run = Polyphemus::Run.create(
          name: workflow_name,
          run_id: workflow_uid,
          config_id: config.config_id,
          version_number: config.version_number,
          created_at: Time.now,
          updated_at: Time.now
        )
        return run
      else
        logger.warn("Failed to submit Argo workflow: #{stderr}")
        raise Etna::Error, stderr
      end
    end

    def get_workflow_status(run)
      cmd = [
        "argo", "get",
        "--output=json",
        run.name,
        "-n", "argo"
      ]

      logger.info("Getting Argo workflow status for #{run.name}")
      logger.info("Argo command: #{cmd.join(" ")}")

      stdout, stderr, status = Open3.capture3(cmd.join(" "))

      return "invalid" unless status.success?

      workflow_data = JSON.parse(stdout)

      return workflow_data["status"]["phase"]
    end

    module_function :submit_workflow, :get_workflow_status
  end
end
