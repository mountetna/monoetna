module Etna
  class Injection
    # Extend into class
    module FromHash
      def from_hash(hash, hash_has_string_keys, rest: nil, key_rest: nil)
        ::Etna::Injection.inject_new(self, hash ,hash_has_string_keys, rest: rest, key_rest: key_rest) do |missing_p|
          raise "required argument '#{missing_p}' of #{self.name} is missing!"
        end
      end
    end

    def self.inject_new(cls, hash, hash_has_string_keys, rest: nil, key_rest: nil, &missing_req_param_cb)
      args, kwds = prep_args(
        cls.method(:initialize), hash, hash_has_string_keys,
        rest: rest, key_rest: key_rest, &missing_req_param_cb
      )
      cls.new(*args, **kwds)
    end

    def self.prep_args(method, hash, hash_has_string_keys, rest: nil, key_rest: nil, &missing_req_param_cb)
      new_k_params = {}
      new_p_params = []

      method.parameters.each do |type, p_key|
        h_key = hash_has_string_keys ? p_key.to_s : p_key
        if type == :rest && rest
          new_p_params.append(*rest)
        elsif type == :keyrest && key_rest
          new_k_params.update(key_rest)
        elsif type == :req || type == :opt
          if hash.include?(h_key)
            new_p_params << hash[h_key]
          elsif type == :req
            if block_given?
              new_p_params << missing_req_param_cb.call(h_key)
            else
              new_p_params << nil
            end
          end
        elsif type == :keyreq || type == :key
          if value.include?(h_key)
            new_k_params[p_key] = hash[h_key]
          elsif type == :keyreq
            if block_given?
              new_k_params[p_key] = missing_req_param_cb.call(h_key)
            else
              new_k_params[p_key] = nil
            end
          end
        end
      end

      [new_k_params, new_p_params]
    end
  end
end