class Magma
  class Rule < Sequel::Model
    one_to_many :identifiers
  end
end
