require 'shellwords'
require 'set'


class Vulcan
  class Snakemake

    class CommandBuilder
      attr_accessor :options, :targets

      # Initialize the CommandBuilder with default options
      def initialize
        @targets = []
        @options = {
          local: false,
          config_path: '',
          profile_path: '',
          workflow_profile_path: '',
          singularity_args: '',
          dry_run: false,
        }
      end

        def build
        # To run any snakemake pipeline, snakemake requires all config values to be present, even if you aren't running
        # that particular job or target. Therefore all snakemake workflows should have a config.yaml with default values.
        #
        # When we pass a custom config each time via --configfile, snakemake overwrites these defaults.
        # BUT in the case the custom config doesn't have all the original params in config.yaml, the remaining
        # default params are still active for the workflow. Therefore we must avoid triggering downstream jobs with these
        # default params.
        #
        # Additionally the Vulcan UI operates in such a way where it just sends params and expects jobs associated to those params to
        # run. Similarly, we only want these jobs to run and not have the default params trigger downstream jobs.
        #
        # Therefore to address both of these concerns we must make sure only the jobs associated with a particular set of
        # params are triggered. However, it is much easier to infer what outputs to trigger and more canonical to
        # actually tell snakemake what outputs or targets to create. This is the entire reason we maintain the target
        # mapping in the db.

        cmd = ['snakemake']
        cmd << "--configfile #{@options[:config_path]}"
        cmd << "--profile #{@options[:profile_path]}"

        # Optional flags and arguments
        cmd << "--workflow-profile #{@options[:workflow_profile_path]}" unless @options[:workflow_profile_path].nil?
        cmd << "--cores 1" if @options[:local]
        unless @options[:singularity_args].nil?
          cmd << "--use-singularity"
          cmd << "--singularity-args '#{Shellwords.escape(@options[:singularity_args])}'"
        end

        if @options[:dry_run]
          cmd << "--dry-run"
        end

        # Targets at the end
        unless targets.empty?
          cmd.concat(targets.map { |target| Shellwords.escape(target) })
        end

        cmd.join(' ')
      end
    end
  end
end
