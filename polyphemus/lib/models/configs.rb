class Polyphemus
  class Config < Sequel::Model

    plugin :timestamps, update_on_create: true, allow_manual_update: true
    plugin :column_encryption do |enc|
      enc.key 0, Polyphemus.instance.config(:etl_secrets_key)
      enc.column :secrets, format: :json
    end

    def self.current
      return self.reverse(:config_id, :version_number).distinct(:config_id).from_self
    end


    def self.next_id
      Polyphemus.instance.db.get { nextval("etl_configs_ids") }
    end

    def as_json(options = nil)
      to_hash.map do |key,value|
        case key
        when :secrets
          next [ key, secrets.map do |skey, secret|
            [ skey, secret.to_s.empty? ? '' : '***' ]
          end.to_h ]
        when :id
          next
        else
          case value
          when Time
            next [ key, value.iso8601 ]
          end

          next [ key, value ]
        end
      end.compact.to_h
    end

    def to_revision
      as_json.slice(:config, :created_at, :version_number, :comment)
    end

  end
end
