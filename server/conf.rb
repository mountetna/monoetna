# conf.rb
# Configuration for Metis.

module Conf

  # Names and order of the fields to be signed.
  SIGNATURE_ITEMS = [

    'directory',
    'hashing_algorithm',
    'start_timestamp',
    'token',
    'original_name', 
    'file_name',
    'file_size',
    'user_email',
    'user_id',
    'project_id',
    'group_id',
    'db_index'
  ]

  # This is the root directory for the data.
  ROOT_DIR = '/data'

  JANUS_ADDR = 'http://janus-dev.ucsf.edu'

  # One day to finish upload.
  UPLOAD_EXPIRE = 86400

  # Warning and error messages that will end up in the log
  WARNS = [

    :PARAMS_NOT_PRESENT,   # 0
    :PARAMS_NOT_CORRECT,   # 1
    :USER_TOKEN_NOT_VALID, # 2
    :NO_PERMS,             # 3
    :INVALID_GROUP,        # 4
    :FILE_EXISTS,          # 5
    :STATUS_EXISTS,        # 6
    :REQUEST_INVALID,      # 7
    :SENT_BLOB_MALFORMED,  # 8
    :USER_NOT_AUTHORIZED,  # 9
    :NO_FILE_EXISTS,       # 10
    :FILE_REMOVE_FAILED    # 11
  ]

  ERRORS = [

    :FILE_SYS_NOT_SYNCED,   # 0
    :MISSING_STATUS,        # 1
    :FILE_NO_STATUS,        # 2
    :TEMP_NO_STATUS,        # 3
    :TEMP_AND_FILE,         # 4
    :STATUS_NO_TEMP_OR_FILE,# 5
    :CONNECTION_ERROR,      # 6
    :JANUS_SERVER_ERROR     # 7
  ]
end