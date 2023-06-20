class Magma
  class Flag < Sequel::Model
  end

  module Flags
    # All hashes nested under this namespace represent flags that can be used in a production.
    # Each has should define a name key representing the name of a flag and each
    # hash should also specify all potential values that can be assigned to a flag.
    GNOMON_MODE = {
      name: 'gnomon_mode',
      identifier: 'identifier',
      pattern: 'pattern',
      none: 'none'
    }
  end
end