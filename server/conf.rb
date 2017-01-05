# conf.rb
# Configuration for Metis.

module Conf

  # This is a temp secret key DO NOT USE THIS IN PRODUCTION!
  SECRET_KEY = 'eOHJ6a48xyOmmaYIM0AnUAOIiJZ0kOZzCxTVtSrhagg'

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
    'user_email',
    'user_id',
    'project_id',
    'redis_index'
  ]

  # This is the root directory for the data.
  ROOT_DIR = '/data'

  # Upload state errors.
  # FILE_NO_STATUS = 0 # A file exsists, but does not have a status entry.
  # TEMP_NO_STATUS = 1 # A temp file exsists, but does not have a status entry.
  # TEMP_AND_FILE  = 3 # Both the temp file and regular file exist.
  # STATUS_NO_TEMP_OR_FILE = 4 # A status exsits but no corresponding file.

  APP_KEY = 'eiE0XQBb33saliJlrrlQJl1UYlUrqD49CHNyGfFozoo'
  JANUS_ADDR = 'http://janus-dev.ucsf.edu'

  # One day to finish upload.
  UPLOAD_EXPIRE = 86400
end