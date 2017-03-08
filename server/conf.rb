module Conf

  JANUS_ADDR = 'http://janus-dev.ucsf.edu'

  # Warning and error messages that will end up in the log
  WARNS = [

    :PARAMS_NOT_PRESENT, # 0
    :NOT_ADMIN_USER      # 1
  ]

  ERRORS = [

    :JANUS_SERVER_ERROR, # 0
    :CONNECTION_ERROR    # 1
  ]
end