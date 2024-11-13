class Vulcan
  class Config < Sequel::Model(:configs)
    many_to_one :workspaces, :key => :workspace_id
    one_to_many :runs
  end
end
