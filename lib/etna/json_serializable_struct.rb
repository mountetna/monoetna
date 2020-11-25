module Etna
  module JsonSerializableStruct
    def self.included(cls)
      cls.instance_eval do
        def self.as_json(v)
          if v.respond_to? :as_json
            return v.as_json
          end

          if v.is_a? Hash
            return v.map { |k, v| [k, as_json(v)] }.to_h
          end

          if v.class.include? Enumerable
            return v.map { |v| as_json(v) }
          end

          v
        end
      end
    end

    def as_json(keep_nils: false)
      inner_json = members.map do |k|
        v = self.class.as_json(send(k))
        [k, v]
      end.to_h

      return inner_json if keep_nils
      inner_json.delete_if { |k, v| v.nil? }
    end

    def to_json
      as_json.to_json
    end
  end
end
