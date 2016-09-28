# conf.rb
# Configuration for Metis.

module Conf

  # This is a temp secret key DO NOT USE THIS IN PRODUCTION!
  SECRET_KEY = '3b43cdb6a85ce66879620a4ff9cf5039'

  # Names and order of the fields to be signed.
  SIGNATURE_ITEMS = [

    'directory',
    'expires',
    'signing_algorithm',
    'hashing_algorithm',
    'start_timestamp',
    'authorization_token',
    'original_name', 
    'file_name',
    'file_size',
    'user_id',
    'group_id'
  ]

  # This is the root directory for the data.
  ROOT_DIR = '/data'

  # Upload state errors.
  # FILE_NO_STATUS = 0 # A file exsists, but does not have a status entry.
  # TEMP_NO_STATUS = 1 # A temp file exsists, but does not have a status entry.
  # TEMP_AND_FILE  = 3 # Both the temp file and regular file exist.
  # STATUS_NO_TEMP_OR_FILE = 4 # A status exsits but no corresponding file.
end