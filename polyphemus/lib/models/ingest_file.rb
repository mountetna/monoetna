class Polyphemus
  class IngestFile < Sequel::Model
    plugin :timestamps, update_on_create: true
  end
end
