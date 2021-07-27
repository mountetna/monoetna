require "etna"

class MetisGitHelper
  METIS_PATH_REGEX =

    def self.config_file_path
      ::File.join(::Dir.home, "etna.yml")
    end

  def initialize(config, environment = :production)
    raise "Make sure you have defined #{self.config_file_path}." unless config && config.has_key?(environment)

    @config = config
    @environment = environment
  end

  def token
    @token = ENV["TOKEN"]

    if !@token
      stderr("No environment variable TOKEN is set. Set your token with `export TOKEN=<your.janus.token>` before using this helper.")
      exit(1)
    end

    @token
  end

  def metis_client
    @metis_client ||= Etna::Clients::Metis.new(
      token: token,
      host: @config[@environment][:metis][:host],
    )
  end

  def run
    while true
      line = readline()
      puts line
      case line
      when "capabilities"
        # stdout("option")
        stdout("list")
        stdout("fetch")
        # stdout("push")
        # when /^option/
        #   #blah
      when /^list/
        list_refs(line)
        # when /^push/
        # push
      when /^fetch/
        fetch_objects
      when ""
        break
      else
        stderr("Unsupported: #{line}")
        exit(1)
      end
    end
  end

  private

  def stdout(message = "\n")
    STDOUT.puts(message)
  end

  def stderr(message = "\n")
    STDERR.puts(message)
  end

  def readline
    STDIN.readline.strip
  end

  def list_refs(line)
    puts line
    refs.each do |ref|
      stdout("#{ref.sha} #{ref.ref}")
    end
    stdout("@#{head_ref.ref} HEAD")
    stdout
  end

  def refs
    [Ref.new("fake-sha", "fake-ref")]
  end

  def head_ref
    # How do we generate a version and content for HEAD?
    Ref.new("fake-sha", "fake-ref")
  end

  class Ref
    def initialize(sha, ref)
      @sha = sha
      @ref = ref
    end

    def sha
      @sha
    end

    def ref
      @ref
    end
  end
end
