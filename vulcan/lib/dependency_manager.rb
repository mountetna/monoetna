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

    def target_image(interpreter, session)
      # Specify the image with sha here
      target = Vulcan.instance.config(:archimedes_interpreters)&.dig(interpreter.to_sym) || Vulcan.instance.config(:archimedes_run_image)

      dependency_name_w_sha(session, target)
    end

    def interpreter(script:)
      first_line = script.split("\n").first

      # Ignore shebang variations for now
      return "node" if first_line == "#!/usr/bin/env node"

      "python"
    end

    def archimedes_run_sha(session)
      dependency_name_w_sha(session, Vulcan.instance.config(:archimedes_run_image))
    end

    private

    # def git_sha

    # end

    def dependency_name_w_sha(session, image_name)
      # If the session's reference figure has dependencies,
      #   we'll try to honor those in the image passed back.
      # Otherwise, default is just what is in Vulcan config.
      base_image = image_name.split(":").first.to_s

      begin
        sha = dependency_sha(session, base_image)

        return "#{base_image}@#{sha}" unless sha.nil?
      end if session.dependencies

      image_name
    end

    def dependency_sha(session, dependency_name)
      return nil unless session.dependencies.keys.include?(dependency_name)

      sha = session.dependencies[dependency_name]
      return nil if sha.empty?

      sha
    end

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
