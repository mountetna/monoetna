class Metis
  class DataBlockMetadata < Sequel::Model
    plugin :timestamps, update_on_create: true
    many_to_one :data_block
  end
end
