class Magma
  class Flag < Sequel::Model
  end

  module Flags
    # All hashes nested under this namespace represent flags that can be used in a production.
    # Each hash should define a 'name' key, that represents the name of the flag and it
    # should also specify all potential values that can be assigned to a flag.
    GNOMON_MODE = {
      name: 'gnomon_mode',
      identifier: 'identifier',
      pattern: 'pattern',
      none: 'none'
    }
  end
end