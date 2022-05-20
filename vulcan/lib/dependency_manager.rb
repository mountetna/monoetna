class Vulcan
  class DependencyManager
    def dependency_shas
      dependencies.map do |dependency|
        [dependency, image_sha(image_name: dependency)]
      end.to_h
    end

    def dependencies
      @dependencies ||= ["vulcan"].concat(
        Vulcan.instance.config(:archimedes_interpreters)&.values&.map { |ai| image_wo_tag(ai) } || []
      ).concat([image_wo_tag(Vulcan.instance.config(:archimedes_run_image))]).uniq.compact
    end

    def tag(dependency = ":production")
      dependency
    end

    def target_image(interpreter)
      Vulcan.instance.config(:archimedes_interpreters)&.dig(interpreter.to_sym) || Vulcan.instance.config(:archimedes_run_image)
    end

    def interpreter(script:)
      first_line = script.split("\n").first

      # Ignore shebang variations for now
      return "node" if first_line == "#!/usr/bin/env node"

      "python"
    end

    def archimedes_run_sha(session)
      # If the session's reference figure has dependencies,
      #   we'll try to honor those in the image passed back.
      # Otherwise, default is just what is in Vulcan config.
      config_image = Vulcan.instance.config(:archimedes_run_image)

      base_image = config_image.split(":").first.to_s

      return "#{base_image}@#{session.dependencies[base_image]}" if session.dependencies && session.dependencies.keys.include?(base_image)

      config_image
    end

    private

    def image_wo_tag(image_name)
      image_name.to_s.split(":").first
    end

    def image_sha(image_name:, postfix: ":production", prefix: "etnaagent/")
      # This assumes that the images are available locally and up to date.
      # Alternative is to ping docker hub or some other registry, but
      #   I'm hesitant to introduce a runtime dependency on a third-party service...
      image_name_full = image_name
      image_name_full = "#{prefix}#{image_name_full}" unless image_name_full.start_with?(prefix)
      image_name_full = "#{image_name_full}#{postfix}" unless image_name_full.end_with?(postfix)

      image_name_wo_sha = "#{image_name}@"
      image_name_wo_sha = "#{prefix}#{image_name_wo_sha}" unless image_name_wo_sha.start_with?(prefix)
      digest = `docker inspect --format '{{index .RepoDigests 0}}' #{image_name_full}`
      digest.gsub(image_name_wo_sha, "").strip
    end
  end
end
