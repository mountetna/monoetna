require 'shellwords'

class Vulcan
  class SnakemakeCommandBuilder
    attr_accessor :options

    def initialize
      @options = {
        run_until: '',
        local: false,
        config_path: '',
        profile_path: '',
        workflow_profile_path: '',
        singularity_args: ''
      }
    end

    def build
      cmd = ['snakemake']
      # Mandatory
      cmd << "--configfile #{@options[:config_path]}"
      cmd << "--profile #{@options[:profile_path]}"
      cmd << "--until #{@options[:run_until]}"
      # Optional
      cmd << "--workflow-profile #{@options[:workflow_profile_path]}" unless @options[:workflow_profile_path].nil?
      cmd << "--cores 1" if @options[:local]
      unless @options[:singularity_args].nil?
        cmd << "--use-singularity"
        cmd << "--singularity-args '#{Shellwords.escape(@options[:singularity_args])}'"
      end
      cmd.join(' ')
    end
  end

end