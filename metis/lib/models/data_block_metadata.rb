class Metis
  class DataBlockMetadata < Sequel::Model(:data_block_metadata)
    plugin :timestamps, update_on_create: true
    many_to_one :data_block
  end
end
