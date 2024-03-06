class Vulcan
  class Workflow < Sequel::Model
    one_to_many :workspaces
  end
end
