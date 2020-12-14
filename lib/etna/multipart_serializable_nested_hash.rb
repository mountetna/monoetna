module Etna
  module MultipartSerializableNestedHash
    def self.included(cls)
      cls.instance_eval do
        def self.encode_multipart_pairs(value, base_key, is_root, &block)
          if value.is_a? Hash
            value.each do |k, v|
              encode_multipart_pairs(v, is_root ? k : "#{base_key}[#{k}]", false, &block)
            end
          elsif value.is_a? Array
            value.each_with_index do |v, i|
              # This is necessary to ensure that arrays of hashes that have hetergenous keys still get parsed correctly
              # Since the only way to indicate a new entry in the array of hashes is by re-using a key that existed in
              # the previous hash.
              if v.is_a? Hash
                encode_multipart_pairs(i, "#{base_key}[][_idx]", false, &block)
              end

              encode_multipart_pairs(v, "#{base_key}[]", false, &block)
            end
          else
            raise "base_key cannot be empty for a scalar value!" if base_key.length == 0

            if value.respond_to?(:read)
              yield [base_key, UploadIO.new(value, 'application/octet-stream'), {filename: 'blob'}]
            else
              yield [base_key, value.to_s]
            end
          end
        end


        def self.encode_multipart_content(value, base_key = '', is_root = true)
          result = []
          self.encode_multipart_pairs(value, base_key, is_root) { |pair| result << pair }
          result
        end
      end
    end

    def encode_multipart_content(base_key = '')
      value = self
      if value.respond_to? :as_json
        value = value.as_json
      end

      self.class.encode_multipart_content(value, base_key)
    end
  end
end
