# Run the Archimedes checker on all files in the workflow/scripts/
#   directory, to make sure they meet the DSL.

describe "Vulcan scripts" do
  def command(filepath)
    [
      "docker",
      "run",
      "--rm",
    ] + docker_run_args(filepath)
  end

  def archimedes_path(filepath)
    "/app/vulcan_scripts/#{::File.basename(filepath)}"
  end

  def scripts_path
    ::File.join(
      ::File.dirname(__FILE__),
      "..",
      "lib",
      "server",
      "workflows",
      "scripts"
    )
  end

  def mount_path
    path = scripts_path

    if (host_dir_map = ENV["HOST_DIR_MAP"])
      container_path, host_dir = host_dir_map.split("=", 2)
      path.sub!(/^#{container_path}/, host_dir)
    end

    path
  end

  def docker_run_args(filepath)
    [
      "-i",
      "-v",
      "/var/run/docker.sock:/var/run/docker.sock:ro",
      "-v",
      "#{Vulcan.instance.config(:archimedes_exec_volume)}:/archimedes-exec",
      "-v",
      "#{mount_path}:/app/vulcan_scripts",
      Vulcan.instance.config(:archimedes_run_image),
      "poetry",
      "run",
      "archimedes-checker",
      archimedes_path(filepath),
    ]
  end

  def check_script(filepath)
    status = nil
    output_str = nil

    cmd = command(filepath)
    Open3.popen2(*cmd) do |input, output, wait_thr|
      input.close
      output_str = output.read
      status = wait_thr.value.exitstatus
    end

    output_str
  end

  it "all pass Archimedes checker" do
    errors = []
    ::Dir["#{scripts_path}/*.py"].each do |path|
      output = check_script(path)
      errors << output unless "" == output
    end

    expect(errors).to eq([])
  end
end
