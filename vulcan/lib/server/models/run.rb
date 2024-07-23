class Vulcan
  class Run < Sequel::Model
    many_to_one :workspace

    def run_config_path
      "#{workspace.path}/#{config_path}"
    end

    def snakemake_log_path
      "#{workspace.path}/#{log_path}"
    end

  end
end
