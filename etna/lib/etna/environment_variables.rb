module Etna
  module EnvironmentVariables
    # a <- b
    def self.deep_merge(a, b)
      if a.is_a?(Hash)
        if b.is_a?(Hash)
          b.keys.each do |b_key|
            a[b_key] = deep_merge(a[b_key], b[b_key])
          end

          return a
        end
      end

      a.nil? ? b : a
    end

    def load_from_env(prefix, root: {}, env: ENV, downcase: true, sep: '__', &path_to_value_mapper)
      env.keys.each do |key|
        next unless key.start_with?(prefix)

        path = key.split(sep, -1)
        if downcase
          path.each(&:downcase!)
        end

        result = path_to_value_mapper.call(path, env[key])
        next unless result

        path, value = result

        if path.empty?
          root = EnvironmentVariables.deep_merge(root, value)
          next
        end

        # path.shift # drop the prefix.
        # path = key.sub(/_FILE/, '').split("__", -1)
        # path.shift # drop the first, just app name

        target = root
        while path.length > 1
          n = path.shift
          target = (target[n] ||= {})
        end

        target[path.last] = deep_merge(target[path.last], value)
      end
    end
  end
end