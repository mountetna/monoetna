class Group < Sequel::Model
  one_to_many :projects
end
