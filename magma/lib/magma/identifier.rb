class Magma
  class Identifier < Sequel::Model
     many_to_one :grammar
     many_to_one :renamed_to, class: self, key: :renamed_to_id
     one_to_many :aliases, class: self, key: :renamed_to_id
  end
end
