class Magma
  class Flag < Sequel::Model
    # update_or_create is natively implemented in a later version of Sequel - I believe it can also be installed as a plugin
    def self.update_or_create(update_hash, value_hash)
      # handles updates, creates and deletes. Deletes when {value: nil}.
      #
      # update_hash: a hash with two keys: {project_name: , flag_name: }
      # value: a hash with one key: {value: }

      existing_flag = self.where(update_hash).first
      if existing_flag
        # We delete flags set as nil and update the flags otherwise
        value_hash[:value].nil? ? existing_flag.delete : existing_flag.update(value_hash)
      else
        # create a new row, make sure the value is not nil
        # our unless statement is for a rare case when a user tries to delete a row that does not exist
        self.create(update_hash.merge(value_hash)) unless value_hash[:value].nil?
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

    def self.is_registered(flag_name, flag_value)
      # flag_hash: a string, flag_value: a string
      # ensures that a flag has the proper name and the values are defined
      constants.each do |constant|
        next unless const_get(constant)[:name] == flag_name && (const_get(constant).has_value?(flag_value) || flag_value.nil?)
        return true
      end
      false
    end

  end
end