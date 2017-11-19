class Upload < Sequel::Model
  one_to_one :file
end
