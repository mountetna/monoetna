require_relative '../../etls/metis/loader'

class Polyphemus
  class MetisLinkerManifest < Polyphemus::WorkflowManifest

    def self.as_json
      {
        name: 'metis',
        schema: Metis::Loader.to_schema,
        runtime_params: {
          commit: {
            type: 'boolean',
            description: 'Commit results to Magma',
            default: false
          },
          updates_only: {
            type: 'boolean',
            description: 'Only link files updated on Metis since the last loader run',
            default: false
          },
          debug: {
            type: 'boolean',
            description: 'Include more information in the run log, e.g. all file paths even when there are many files',
            default: false
          }
        },
        workflow_path: '/app/workflows/argo/metis_linker/workflow.yaml'
      }
    end

    def validate
    end

    private

    def array_or_string_param(param, allowed_values=["all"])
      request_params[param].is_a?(Array) || allowed_values.include?(request_params[param])
    end

  end
end
