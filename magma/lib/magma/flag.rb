class Magma
  class Flag < Sequel::Model
    # update_or_create is natively implemented in a later version of Sequel - I believe it can also be installed as a plugin
    def self.update_or_create(project_name, flag_list)
      flag_list.each do |flag|
        flag.each do |flag_name, value|
          existing_flag = self.where(project_name: project_name, flag_name: flag_name.to_s).first
          if existing_flag
            existing_flag.update(value: value)
          else
            self.create(project_name: project_name, flag_name: flag_name.to_s, value: value)
          end
        end
      end
    end
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