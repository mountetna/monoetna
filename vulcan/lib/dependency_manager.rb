class Vulcan
  class DependencyManager
    def initialize
      @existing_images_cache = {}
    end

    def dependency_shas
      dependencies.map do |dependency|
        [dependency, image_sha(image_name: dependency)]
      end.concat([git_sha]).to_h
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

      return target if ignore_dependencies?

      tag_or_latest_image(dependency_name_w_sha(session, target))
    end

    def interpreter(script:)
      first_line = script.split("\n").first

      # Ignore shebang variations for now
      return "node" if first_line == "#!/usr/bin/env node"
      # Require a library call in first line as hint for R scripts.
      return "r" if first_line.include? "library(" 
      return "r" if first_line.include? "::load_packages"

      "python"
    end

    def archimedes_run_sha(session)
      configured_archimedes = Vulcan.instance.config(:archimedes_run_image)

      return configured_archimedes if ignore_dependencies?

      tag_or_latest_image(
        dependency_name_w_sha(session, configured_archimedes)
      )
    end

    def ignore_dependencies?
      @ignore_dependencies ||= begin
        if [:development, :test].include?(Vulcan.instance.environment)
          return !!Vulcan.instance.config(:ignore_dependencies)
        end

        false
      end
    end

    private

    def git_sha
      ['monoetna_git_sha', ENV['MONOETNA_SHA'] || '']
    end

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

    def tag_or_latest_image(image_name)
      # If the image is specified with a SHA (@sha:<hash>), then
      #   check that the hash exists in the registry. If not,
      #   then return the image name with :latest tag, as the
      #   fallback.
      image_exists?(image_name) ? image_name : latest_image(image_name)
    end

    def latest_image(image_name)
      # Convert etnaagent/image@sha256:<hash> or etnaagent/image:<tag>
      # to etnaagent/image:latest
      "#{image_name.split("@").first.split(":").first}:latest"
    end

    def image_exists?(image_name)
      key_check_result = @existing_images_cache[image_name]
      return key_check_result unless key_check_result.nil?

      # Run a docker command to check that the sha or image exists.
      exists = `docker manifest inspect #{image_name} > /dev/null`

      @existing_images_cache[image_name] = 0 == $?.exitstatus ? true : false
    end
  end
end
