# conf.rb
# Configuration for Metis.

module Conf

  # This is a temp secret key DO NOT USE THIS IN PRODUCTION!
  SECRET_KEY = '<%= metis_secret_key %>?'

  # Names and order of the fields to be signed.
  SIGNATURE_ITEMS = [ <%= metis_signature_items%> ]

  # This is the root directory for the data.
  ROOT_DIR = '/data'

  # Upload state errors.
  # FILE_NO_STATUS = 0 # A file exsists, but does not have a status entry.
  # TEMP_NO_STATUS = 1 # A temp file exsists, but does not have a status entry.
  # TEMP_AND_FILE  = 3 # Both the temp file and regular file exist.
  # STATUS_NO_TEMP_OR_FILE = 4 # A status exsits but no corresponding file.

  APP_KEY = '<%= metis_app_key %>'
  JANUS_ADDR = '<%= janus_url %>'
end