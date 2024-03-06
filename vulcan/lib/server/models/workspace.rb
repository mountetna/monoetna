class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow
    one_to_many :runs
  end
end
