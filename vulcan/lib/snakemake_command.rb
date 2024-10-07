require 'shellwords'
require 'set'


class Vulcan
  class Snakemake

    class CommandBuilder
      attr_accessor :target_meta, :options

      # Initialize the CommandBuilder with default options
      def initialize
        @target_meta = {
          mapping: {},
          provided_params: [],
          available_files: [],
        }
        @options = {
          local: false,
          config_path: '',
          profile_path: '',
          workflow_profile_path: '',
          singularity_args: '',
          dry_run: false,
        }
      end

      def find_buildable_targets(target_mapping, provided_params_set, available_files_set)
        buildable = Set.new
        target_mapping.each do |target, requirements|
          # Check if all input files are available
          inputs_satisfied = requirements["inputs"].to_set.subset?(available_files_set)

          # Check if all required params are provided
          params_satisfied = requirements["params"].to_set.subset?(provided_params_set)

          # If both conditions are satisfied, add the target to buildable
          if inputs_satisfied && params_satisfied
            buildable.add(target)
          end
        end

        buildable
      end


      def determine_buildable_targets
        all_targets = Set.new
        all_files = @target_meta[:available_files].to_set
        loop do
          new_targets = find_buildable_targets(@target_meta[:mapping], @target_meta[:provided_params].to_set, all_files)
          # Are there any new targets?
          break if all_targets == new_targets
          all_targets.merge(new_targets)
          all_files.merge(new_targets) # we add these targets to files that can be generated
        end
        all_targets.to_a
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
        targets = determine_buildable_targets
        unless targets.empty?
          cmd.concat(targets.map { |target| Shellwords.escape(target) })
        end

        cmd.join(' ')
      end
    end
  end
end
