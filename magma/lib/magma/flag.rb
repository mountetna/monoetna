class Magma
  class Flag < Sequel::Model
  end

  module Flags
    # All modules nested under this namespace represent flags that can be used in a production.
    # Each module should define a NAME constant representing the name of a flag and each
    # module should also specify all potential values that can be assigned to a flag.
    #
    GNOMON_MODE = {
      name: 'gnomon_mode',
      identifier: 'identifier',
      pattern: 'pattern',
      none: 'none'
    }
  end
end