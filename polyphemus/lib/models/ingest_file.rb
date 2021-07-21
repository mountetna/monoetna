class Polyphemus
  class IngestFile < Sequel::Model
    plugin :timestamps, update_on_create: true

    def to_hash
      {
        name: name,
        host: host,
      }
    end
  end
end
