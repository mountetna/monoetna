# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()

Metis.add_route('GET', '/', 'ClientController#index')
Metis.add_route('GET', '/admin', 'ClientController#admin')
Metis.add_route('GET', '/user', 'ClientController#user')

Metis.add_route('POST', '/login', 'LogController#log_in')
Metis.add_route('POST', '/logout', 'LogController#log_out')
Metis.add_route('POST', '/check', 'LogController#check_log')

Metis.add_route('POST', '/upload-authorize', 'UploadController#authorize_upload')
Metis.add_route('POST', '/upload-start', 'UploadController#start_upload')
Metis.add_route('POST', '/upload-blob', 'UploadController#upload_blob')
Metis.add_route('POST', '/upload-pause', 'UploadController#pause_upload')
Metis.add_route('POST', '/file-remove', 'UploadController#remove_file')

Metis.add_route('POST', '/retrieve-files', 'Controller#retrieve_files')

# These items will delegate to Janus. i.e. we basically make pass through
# requests to janus.
Metis.add_route('POST', '/get-users', 'AdminController#get_users')
Metis.add_route('POST', '/get-projects', 'AdminController#get_projects')
Metis.add_route('POST', '/get-permissions', 'AdminController#get_permissions')
Metis.add_route('POST', '/save-permission', 'AdminController#save_permission')
Metis.add_route('POST', '/upload-permissions', 'AdminController#upload_permissions')
