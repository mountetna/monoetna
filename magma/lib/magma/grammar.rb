class Magma
  class Grammar < Sequel::Model
    one_to_many :identifiers
  end
end
