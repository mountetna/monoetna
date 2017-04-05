# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()

Metis.add_route('GET', '/', 'ClientController#index')
Metis.add_route('GET', '/user', 'ClientController#user')

Metis.add_route('POST', '/login', 'UserLogController#log_in')
Metis.add_route('POST', '/logout', 'UserLogController#log_out')
Metis.add_route('POST', '/check', 'UserLogController#check_log')

Metis.add_route('POST', '/upload-authorize', 'UploadController#authorize_upload')
Metis.add_route('POST', '/upload-start', 'UploadController#start_upload')
Metis.add_route('POST', '/upload-blob', 'UploadController#upload_blob')
Metis.add_route('POST', '/upload-pause', 'UploadController#pause_upload')
Metis.add_route('POST', '/upload-cancel', 'UploadController#cancel_upload')
Metis.add_route('POST', '/remove-file', 'UploadController#remove_file')
Metis.add_route('POST', '/remove-failed', 'UploadController#remove_failed')

Metis.add_route('POST', '/retrieve-files', 'MainController#retrieve_files')