class Vulcan
  class Run < Sequel::Model
    many_to_one :workspace
  end
end
