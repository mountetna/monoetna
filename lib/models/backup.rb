class Metis
  class Backup < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_many :files
  end
end
