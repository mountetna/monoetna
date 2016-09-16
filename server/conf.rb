# Configuration for Metis

module Conf

  # This is a temp secret key DO NOT USE THIS IN PRODUCTION
  SECRET_KEY = '3b43cdb6a85ce66879620a4ff9cf5039'

  # Names and order of the fields to be signed
  SIGNATURE_ITEMS = [

    'directory', 
    'expires',
    'algorithm',
    'timestamp',
    'type',
    'user_email',
    'auth_token',
    'file_name',
    'file_size'
  ]

  # This is the root directory for the data
  ROOT_DIR = '/data'
end