# conf.rb
# Configuration for Metis.

module Conf

  # Names and order of the fields to be signed.
  SIGNATURE_ITEMS = [

    'start_timestamp',
    'token',
    'file_name',
    'original_name',
    'file_size',
    'user_email',
    'group_name',
    'project_name',
    'hashing_algorithm'
  ]

  FILE_VALIDATION_ITEMS = {

    'original_name'=> String,
    'file_name'=> String,
    'file_size'=> Integer,
    'group_name'=> String,
    'project_name'=> String,
    'start_timestamp'=> Integer,
    'user_email'=> String,
    'hashing_algorithm'=> String
  }

  UPLOAD_VALIDATION_ITEMS = {

    'current_byte_position'=> Integer,
    'current_blob_size'=> Integer,
    'next_blob_size'=> Integer,
    'next_blob_hash'=> String,
    'status'=> String
  }

  # This is the root directory for the data.
  ROOT_DIR = '/data'

  JANUS_ADDR = 'https://janus-stage.ucsf.edu'

  # One day to finish upload.
  UPLOAD_EXPIRE = 86400

  # Warning and error messages that will end up in the log
  WARNS = [

    :PARAMS_NOT_PRESENT,   # 0
    :JANUS_CONN_ERROR,     # 1
    :JANUS_CONN_FAIL,      # 2
    :USER_TOKEN_NOT_VALID, # 3
    :FILE_EXISTS,          # 4
    :METADATA_EXISTS,      # 5
    :MISSING_FILE,         # 6
    :MISSING_METADATA,     # 7
    :REQUEST_INVALID,      # 8
    :SENT_BLOB_MALFORMED   # 9
  ]

  ERRORS = [

    :HTTPS_REQ_FAILED,     # 0
    :JANUS_CONN_ERROR,     # 1
    :FILE_SYS_NOT_SYNCED   # 2
  ]
end